#!/usr/bin/python
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

from __future__ import print_function, division, absolute_import
__version__ = 1.1
__author__ = 'Tim Booth, Edinburgh Genomics <tim.booth@ed.ac.uk>'

import os, sys, re
import struct
import gzip

#For handling byte strings I need have alternative code for Py2 vs. Py3
#so set a global flag.  Running this string check within the loop adds a
#noticeable slow-down.
PY3 = (sys.version >= '3')

#Callers may find these constants useful when dealing with results.
SEQUENCE  = 0
QUAL_FLAG = 1

class BCLReader(object):

    def __init__(self, location="."):
        """Creates a BCLReader instance that reads from a single run.
           location: The top level data directory for the run.
           This should be the one that contains the Data directory and the
           RunInfo.xml file.
        """
        #Just check that we can read the expected files at this
        #location.
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

        #return nuc_string, flag
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

        #Find the file prefix I need to be looking at.  Could infer it
        #from the lane number but instead I'll do it by looking for the
        #matching .filter file
        self.prefix = None
        self.data_dir = data_dir
        self.tile = tile

        data_dir_listing = os.listdir(data_dir)
        for filt in data_dir_listing:
            amatch = re.match('(.+_%s).filter' % tile, filt)
            if amatch:
                self.prefix = amatch.group(1)
                self.filter_file = os.path.join(data_dir, filt)
                break

        if not self.prefix:
            raise RuntimeError("Cannot find a .filter file for tile %s" % tile)

        #The CBCL profix is in the format "L00{lane}_{surface}", where the prefix
        #should have the same name as the data_dir, and the surface is the first
        #digit of the tile name.
        self.cbcl_prefix = "%s_%s" % (os.path.basename(data_dir), str(tile)[0])

        #Also work out the number of cycle folders.  Should be 308
        #for HiSeq 4000
        cycle_dirs = [ f for f in data_dir_listing if re.match('C\d+.1$', f) ]
        self.num_cycles = len(cycle_dirs)

        #And also the number of clusters, which should be 4309650
        #for HiSeq 4000.  We need to snag this from the .filter file
        with open(self.filter_file, 'rb') as filt_fh:

            filt_header = struct.unpack('<III', filt_fh.read(12))
            #The first word should be 0, and the version byte is 3, at least
            #on my test files.
            assert tuple(filt_header[0:2]) == (0, 3)
            self.num_clusters = filt_header[2]

    def get_seqs(self, cluster_indices, start=0, end=None):
        """Collects the sequences specified by indices as a hash of pairs of seq+flag.  Ie.
                result = { idx1: ( 'ATCG...', True ), idx2: ('NNNG...', False), ... }
           indices: An iterable that yields integers.  The order is unimportant.
           start: starting base.
           end: ending base.  Note this is as in Python's range(start,end), so
                start=2 and end=10 will skip the first two bases and yield the
                next 8.
        """
        #To build the sequence we have to loop over all the .bcl.gz files in the
        #cycle folders.  These are all named C[num].1 where num is 1-308 (unpadded).
        #The first base of the read will not be in C1.1 because there is an adapter (really??),
        #and also for paired-end reads it doesn't make sense to read all bases as this
        #will run though onto the other end.  However, we can worry about this later.
        if end is None:
            end = self.num_cycles

        #We could use NumPy
        #for more efficient storage but it seems overkill as the number of reads
        #being extracted is comparatively small.
        #This also ensures that all the indices are ints
        seq_collector =  { int(idx) : ['N'] * (end - start) for idx in cluster_indices }
        flag_collector = { int(idx) : None for idx in cluster_indices }
        sorted_keys = sorted(seq_collector.keys())

        #Fail fast if a key is out of range
        if sorted_keys[-1] >= self.num_clusters:
            raise IndexError("Requested cluster %i is out of range.  Highest on this tile is %i." %
                             (sorted_keys[-1], self.num_clusters-1) )

        #And just to be sure, no key should be negative
        if sorted_keys[0] < 0:
            raise IndexError("Requested cluster %i is a negative number." % sorted_keys[0])

        for cycle in range(start, end):
            cycle_dir = os.path.join(self.data_dir, 'C%i.1' % (cycle + 1))

            # Now are we looking at .bcl.gz files or .cbcl files??
            cycle_file = os.path.join(cycle_dir, '%s.bcl.gz' % self.prefix)
            cbcl_file  = os.path.join(cycle_dir, '%s.cbcl' % self.cbcl_prefix)

            try:
                with gzip.open(cycle_file, 'rb') as fh:
                    self._get_seqs_from_fh(fh, seq_collector, sorted_keys)
            except FileNotFoundError:
                # Try the cbcl file. If this fails allow the stack trace which will report both
                # errors.
                # TODO - develop this code. See cbcl_read.py

                with gzip.open(data_handle, 'rb') as fh:
                    self._get_seqs_from_fh(fh, seq_collector, sorted_keys)


        #Now get the accept/reject flag from the .filter file
        #This must exist as we opened it earlier when reading self.num_clusters
        with open(self.filter_file, 'rb') as filt_fh:

            filt_header = filt_fh.read(12)
            assert tuple(struct.unpack('<III', filt_header)) == (0, 3, self.num_clusters)

            for idx in sorted_keys:
                try:
                    filt_fh.seek(idx + 12)
                except:
                    sys.stderr.write("Error seeking ahead at location: %s + 12"%idx)
                    raise
                filt_byte, = struct.unpack('B', filt_fh.read(1))

                flag_collector[idx] = bool( filt_byte & 0b00000001 )


        #Remap the arrays into strings
        # return dict( idx : (nuc_string, flag) )
        return { idx : ( ''.join(seq), flag_collector[idx] ) for idx, seq in seq_collector.items() }

def _get_seqs_from_fh(self, fh, seq_collector, sorted_keys):
    """ Reads from the fh, which is presumably a gzip stream handle, and
        adds the specified seqs to the seq_collector.
        This is intended for internal use only.
        And obviously it can only be called once per fh.
    """
    bcl_header = fh.read(4)

    #The BCL header should be a fixed length depending on the machine type.
    #This assertion checks that it is at least consistent with the filter
    #file for this tile.
    assert struct.unpack('<I', bcl_header)[0] == self.num_clusters

    #I envisaged a a cunning system where we would seek through the file,
    #just reading the chunks we wanted.  Turns out for more than, say,
    #10 reads, it's faster just to slurp the thing.  For over 10000 it's
    #considerably faster!
    if len(sorted_keys) > 10:
        slurped_file = fh.read()

        for idx in sorted_keys:
            if PY3:
                base_byte = slurped_file[idx]
            else:
                base_byte = ord(slurped_file[idx])

            #base = 'N'
            #qual = 0
            if base_byte:
                #The two lowest bits give us the base call
                base = ('A', 'C', 'G', 'T')[base_byte & 0b00000011]

                #And the high bits give us the quality, but we're not using
                #it here.
                #qual = base_byte >> 2

                #Collect the base
                seq_collector[idx][cycle - start] = base
    else:
        for idx in sorted_keys:
            fh.seek(idx + 4)
            #Is reading bytes 1 at a time slow?  I'd imagine that internal
            #cacheing negates any need for chunked reads at this level.
            if PY3:
                base_byte, = fh.read(1)
            else:
                base_byte, = struct.unpack('B', fh.read(1))

            #Copy-paste-ahoy!
            if base_byte:
                base = ('A', 'C', 'G', 'T')[base_byte & 0b00000011]
                seq_collector[idx][cycle - start] = base


