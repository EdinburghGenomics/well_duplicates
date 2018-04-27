#!/usr/bin/env python3
"""
A module to grab sequence reads direct from the .bcl and .filter files
outputted by Illumina.  The motivation is that for some QC tasks we want
to grab a small subsample of reads, and getting these from the FASTQ is
very inefficient, especially once everything is demultiplexed and zipped.

We already have a C implementation of this in the bcl2fastq source and
a Java implementation in Picard Tools, but the world needs Python.

We assume not only the BCL format but also the standard directory
layout as specified in
https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf

This module can take advantage of the fact that the BCL files have a fixed record length
and use seek() to jump to the location where the position of interest is stored.
Unfortunately for GZipped files this does not give much of an advantage, as the file
must be decompressed internally to perform the seek().  For reading several sequences
at once a simple load to memory turns out to be faster, so the reader will do
that.

For max efficiency you should call get_seqs() just once per tile with
all the locations you want to extract.

Synopsis:

   proj = BCLReader("/your/project/dir")
   tile = proj.get_tile(1, 1101)

   all_seqs = tile.get_seqs([70657,70658,70659], start=20, end=40)
   seq1, flag1 = all_seqs[70567]
   seq2, flag2 = all_seqs[70568]
   seq3, flag3 = all_seqs[70569]

   how_many_valid = sum([flag1,flag2,flag3])

On 24th Oct 2017:
We'd also like this module to be able to read .cbcl files, which are concatenated BCL files
(aka. indexed gzip files). Reading these efficiently might require changing the API a little.

"""

__version__ = 1.2
__author__ = 'Tim Booth, Edinburgh Genomics <tim.booth@ed.ac.uk>'

import os, sys, re
import struct
import gzip

# This now works only in Python3 - byte semantics are totally different
assert sys.version >= '3'

# Callers may find these constants useful when dealing with results.
SEQUENCE  = 0
QUAL_FLAG = 1

class BCLReader(object):

    def __init__(self, location="."):
        """Creates a BCLReader instance that reads from a single run.
           location: The top level data directory for the run.
           This should be the one that contains the Data directory and the
           RunInfo.xml file.
        """
        # Just check that we can read the expected files at this
        # location.
        basecalls_dir = os.listdir( os.path.join(location, "Data", "Intensities", "BaseCalls") )
        self.lanes = [ d for d in basecalls_dir if re.match('L\d\d\d$', d) ]

        self.location = location


    def get_seq(self, lane, tile, cluster_index, start=0, end=None):
        """Fetches a single sequence from a specified tile.
           lane: lane number (see get_tile)
           tile: tile number (see get_tile)
           cluster_index: cluster number counting from 0.  To determine
             this from the standard position co-ordinates you need to
             parse the .locs file for the run.
           ** This is going to be very inefficient for fetching multiple sequences. **
           Use get_tile(...).get_seqs([0,2,4,6]) instead.
        """
        tile = self.get_tile(lane, tile)

        result = tile.get_seqs([cluster_index], start, end)

        # return nuc_string, flag
        return result[cluster_index]

    def get_tile(self, lane, tile, in_memory=False):
        """Opens a tile for reading.  You can then call get_seqs() to actually
           fetch the data.
           Lane and tile should be specified as per the Illumina file structure,
           so lanes are 1 to 8 and tiles are eg. [12][12]{01-28} (for HiSeq 4000).
        """
        lane_dir = str(lane)
        if lane_dir not in self.lanes:
            lane_dir = 'L%03d' % int(lane_dir)

        data_dir = os.path.join(self.location, "Data", "Intensities", "BaseCalls", lane_dir)

        if in_memory:
            raise RuntimeError("Preloading into memory not implemented yet")

        return Tile(data_dir, tile)


class Tile(object):

    def __init__(self, data_dir, tile):
        """Fetches sequences from a single tile.
           You would not normally instantiate these directly.  Create a
           BCLReader and call get_tile() instead.
        """

        # Find the file prefix I need to be looking at.  Could infer it
        # from the lane number but instead I'll do it by looking for the
        # matching .filter file
        self.bcl_filename = None
        self.data_dir = data_dir
        self.tile = tile

        data_dir_listing = os.listdir(data_dir)
        for filt in data_dir_listing:
            amatch = re.match('(.+_%s).filter' % tile, filt)
            if amatch:
                self.bcl_filename = amatch.group(1) + '.bcl.gz'
                self.filter_file = os.path.join(data_dir, filt)
                break

        if not self.bcl_filename:
            raise RuntimeError("Cannot find a .filter file for tile %s" % tile)

        # The CBCL profix is in the format "L00{lane}_{surface}", where the prefix
        # should have the same name as the data_dir, and the surface is the first
        # digit of the tile name.
        self.cbcl_filename = "%s_%s.cbcl" % (os.path.basename(data_dir), str(tile)[0])

        # Also work out the number of cycle folders.  Should be 308
        # for HiSeq 4000
        cycle_dirs = [ f for f in data_dir_listing if re.match('C\d+.1$', f) ]
        self.num_cycles = len(cycle_dirs)

        # And also the number of clusters, which should be 4309650
        # for HiSeq 4000.  We need to snag this from the top of the .filter file
        with open(self.filter_file, 'rb') as filt_fh:

            filt_header = struct.unpack('<III', filt_fh.read(12))
            # The first word should be 0, and the version byte is 3, at least
            # on my test files.
            assert tuple(filt_header[0:2]) == (0, 3)
            self.num_clusters = filt_header[2]

        # For now, don't read the rest of the .filter file
        self.filter_offsets = None
        self.passing_wells = None

    def get_seqs(self, cluster_indices, start=0, end=None):
        """Collects the sequences specified by indices as a hash of pairs of seq+flag.  Ie.
                result = { idx1: ( 'ATCG...', True ), idx2: ('NNNG...', False), ... }
           indices: An iterable that yields integers.  The order is unimportant.
           start: starting base.
           end: ending base.  Note this is as in Python's range(start,end), so
                start=2 and end=10 will skip the first two bases and yield the
                next 8.
        """
        # To build the sequence we have to loop over all the .bcl.gz files for the selected tile
        # in the cycle folders.  These are all named C[num].1 where num is 1-308 (unpadded).
        # The first proper base of the read will be in C1.1 as the read is set to start after the
        # adapter, but there may still be an internal barcode (eg. for RAD) so high diversity
        # is not assured.
        # To understand the meaning of all the cycles and what reads they correspond to we do
        # need to look at the run settings.
        if end is None:
            end = self.num_cycles

        # We could use NumPy
        # for more efficient storage but it seems overkill as the number of reads
        # being extracted is comparatively small.
        # This also ensures that all the indices are ints
        seq_collector =  { int(idx) : ['N'] * (end - start) for idx in cluster_indices }
        flag_collector = { int(idx) : None for idx in cluster_indices }
        sorted_keys = sorted(seq_collector.keys())

        # Fail fast if a key is out of range
        if sorted_keys[-1] >= self.num_clusters:
            raise IndexError("Requested cluster %i is out of range.  Highest on this tile is %i." %
                             (sorted_keys[-1], self.num_clusters-1) )

        # And just to be sure, no key should be negative
        if sorted_keys[0] < 0:
            raise IndexError("Requested cluster %i is a negative number." % sorted_keys[0])

        # Get the accept/reject flag from the .filter file
        fo = self._get_filter_offsets()
        for idx in sorted_keys:
            flag_collector[idx] = (fo[idx] != -1)

        # Now the actual basecalls
        for cycle in range(start, end):
            cycle_dir = os.path.join(self.data_dir, 'C%i.1' % (cycle + 1))

            # Now are we looking at .bcl.gz files or NovaSeq .cbcl files??
            cycle_file = os.path.join(cycle_dir, self.bcl_filename)
            cbcl_file  = os.path.join(cycle_dir, self.cbcl_filename)

            try:
                with gzip.open(cycle_file, 'rb') as bcl_fh:
                    self._get_seqs_from_bcl(bcl_fh, cycle - start, sorted_keys, seq_collector)
            except FileNotFoundError:
                # Try the cbcl file. If this fails allow the stack trace which will report both
                # missing files.
                # Note that this does result in opening the same CBCL file again and again
                # for each tile, but each chunk is only unzipped once.
                with open(cbcl_file, 'rb') as fh:
                    self._get_seqs_from_cbcl(fh, cycle - start, sorted_keys, seq_collector)

        # Remap the arrays into strings
        #  return dict( idx : (nuc_string, flag) )
        return { idx : ( ''.join(seq), flag_collector[idx] ) for idx, seq in seq_collector.items() }

    def _get_filter_offsets(self):
        """ Load the filter file, and convert it to a series of offsets. The actual
            offsets are only used when reading excluded CBCL files but the -1 entries
            will indicate the bad wells.
            The file must exist as we opened it earlier when reading self.num_clusters
        """
        # Lazy load
        if self.filter_offsets:
            return self.filter_offsets

        with open(self.filter_file, 'rb') as filt_fh:

            filt_header = filt_fh.read(12)
            # We already saw this!
            assert tuple(struct.unpack('<III', filt_header)) == (0, 3, self.num_clusters)

            # Slurp the whole thing - file length should match the num_clusters value
            # we already know.
            filt_bytes = struct.unpack('<{}B'.format(self.num_clusters), filt_fh.read())

        # Make the table of offsets
        filt_offsets = [-1] * self.num_clusters
        offset = 0
        for n, flag in enumerate(filt_bytes):
            if flag & 0b00000001:
                filt_offsets[n] = offset
                offset += 1

        self.filter_offsets = filt_offsets
        self.passing_wells = offset

        return self.filter_offsets

    def _get_seqs_from_cbcl(self, fh, cycle_idx, sorted_keys, seq_collector):
        """ Reads from the fh to find the appropriate BCL block and then unpacks
            it to extract the basecalls. Deals with excluded/unexcluded flag, requesting
            the filter_offsets as necessary.
            See cbcl_read.py for a more comprehensive version of CBCL reading code.
        """
        # Assume that fh is positioned at the start and read the header...
        # First 12 bytes are fixed fields.
        header_bytes = fh.read(12)
        h_version, h_size, h_basebits, h_qbits, h_bins = struct.unpack('<HIBBI', header_bytes)

        assert h_version == 1
        assert h_size > 32  #Should actually be 5681 for all the current CBCL files
        assert h_basebits == 2
        assert h_qbits == 2 #6 is valid but we don't support it!
        assert h_bins == 4  #implied if h_qbits is 2

        # We don't care about the quality binning info but we do need the tile count,
        # even though we're pretty sure it will always be 352.
        qbin_and_tc_bytes = fh.read((h_bins * 4 * 2) + 4)
        tile_count, = struct.unpack('<I', qbin_and_tc_bytes[-4:])

        # Now I can get all the file offsets which is what I really wanted.
        # Plus the excluded_flag which is the final byte.
        all_offset_bytes = fh.read( tile_count * 16 + 1 )
        excluded_flag = bool(all_offset_bytes[-1])
        t_bcl_offset = h_size

        # A loop is necessary as I have to tot up the csize values to get the seek offset
        tile_as_int = int(self.tile)
        for t in range(tile_count):
            t_number, t_clusters, t_usize, t_csize = struct.unpack('<IIII', all_offset_bytes[t*16:(t+1)*16])

            if t_number == tile_as_int:
                # Found it!
                break

            t_bcl_offset += t_csize

        # We did find it, right?
        assert t_number == tile_as_int

        # Go to the start of the block of tile data and slurp it all (even if
        # I only want 1 or 2 bases - seems pointless to try and optimise the
        # bases < 10 case)
        fh.seek(t_bcl_offset)
        zipdata = gzip.GzipFile(fileobj=fh, mode='rb').read(t_usize)

        if excluded_flag:
            excluded_offsets = self._get_filter_offsets()

        for welln in sorted_keys:

            wellidx = welln
            if excluded_flag:
                wellidx = excluded_offsets[welln]

            if wellidx == -1:
                # Leave the read as an N
                continue

            if wellidx % 2:
                # Take the high bits
                base_byte = zipdata[wellidx//2] >> 4
            else:
                # Take the low bits
                base_byte = zipdata[wellidx//2] & 0b00001111

            # Finally it's the same as for old BCL.
            if base_byte:
                seq_collector[welln][cycle_idx] = ('A', 'C', 'G', 'T')[base_byte & 0b00000011]

    def _get_seqs_from_bcl(self, fh, cycle_idx, sorted_keys, seq_collector):
        """ Reads from the fh, which is presumably a gzip stream handle, and
            adds the specified seqs to the seq_collector.
            This is intended for internal use only.
            And obviously it can only be called once per fh.
        """
        bcl_header = fh.read(4)

        # The BCL header should be a fixed length depending on the machine type.
        # This assertion checks that it is at least consistent with the filter
        # file for this tile.
        assert struct.unpack('<I', bcl_header)[0] == self.num_clusters

        # I envisaged a a cunning system where we would seek through the file,
        # just reading the chunks we wanted.  Turns out for more than, say,
        # 10 reads, it's faster just to slurp the thing.  For over 10000 it's
        # considerably faster!
        if len(sorted_keys) > 10:
            slurped_file = fh.read()

            for idx in sorted_keys:
                base_byte = slurped_file[idx]

                # base = 'N'
                # qual = 0
                if base_byte:
                    # The two lowest bits give us the base call
                    base = ('A', 'C', 'G', 'T')[base_byte & 0b00000011]

                    # And the high bits give us the quality, but we're not using
                    # it here, other than the above test which catches no-calls.
                    # qual = base_byte >> 2

                    #Collect the base
                    seq_collector[idx][cycle_idx] = base
        else:
            for idx in sorted_keys:
                fh.seek(idx + 4)
                # Is reading bytes 1 at a time slow?  I'd imagine that internal
                # cacheing negates any need for chunked reads at this level.
                base_byte, = fh.read(1)

                # Copy-paste-ahoy!
                if base_byte:
                    base = ('A', 'C', 'G', 'T')[base_byte & 0b00000011]
                    seq_collector[idx][cycle_idx] = base


