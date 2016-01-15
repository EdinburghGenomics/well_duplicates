#!/usr/bin/python

from __future__ import print_function, division, absolute_import

"""
A module to grab sequence reads direct from the .bcl and .filter files
outputted by Illumina.  The motivation is that for some QC tasks we want
to grab a small subsample of reads, and getting these from the FASTQ is
very inefficient, especially once they are demultiplexed.

We assume not only the BCL format but also the standard directory
layout as specified in
https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf

This module takes advantage of the fact that the BCL files have a fixed record length
and uses seek() to jump to the location where the position of interest is stored.
Unfortunately for GZipped files this does not give much of an advantage, as the file
must be decompressed internally to perform the seek().
Therefore for max efficiency one should call get_seqs() just once with all the locations
you want to extract.
"""

import os, sys, re
import struct
import gzip

class BCLReader(object):

    def __init__(self, location="."):
        #Just check that we can read the expected files at this
        #location.
        basecalls_dir = os.listdir( os.path.join(location, "Data", "Intensities", "BaseCalls") )
        self.lanes = [ d for d in basecalls_dir if re.match('L\d\d\d$', d) ]

        self.location = location


    def get_seq_by_loc(self, lane, tile, index, start=0, end=None):
        """This is going to be very inefficient for fetching multiple sequences.
           Use get_tile(...).get_seqs([0,2,4,6]) instead.
        """
        tile = self.get_tile(lane, tile)

        result = tile.get_seqs([index], start, end)

        #return nuc_string, flag
        return result[index]

    def get_tile(self, lane, tile, in_memory=False):
        """Opens a tile for reading.  You can now use get_seqs() to actually
           fetch the data.
           Lane and tile should be specified as per the Illumina file structure,
           so lanes are 1-8 and tiles are eg. [12][12]{01-28} (for HiSeq 4000).
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

        #Find the file prefix I need to be looking at.  Could infer it
        #from the lane number but instead I'll do it by looking for the
        #.filter file
        self.prefix = None
        data_dir_listing = os.listdir(data_dir)
        for filt in data_dir_listing:
            amatch = re.match('(.+_%s).filter' % tile, filt)
            if amatch:
                self.prefix = amatch.group(1)
                self.filterfile = filt
                break

        if not self.prefix:
            raise RuntimeError("Cannot find a .filter file for tile %s" % tile)

        self.data_dir = data_dir
        self.tile = tile

        #Also work out the number of cycle folders.  Should be 308
        cycle_dirs = [ f for f in data_dir_listing if re.match('C\d+.1$', f) ]
        self.num_cycles = len(cycle_dirs)

        #TODO - try opening all 307+1 filehandles in advance.  Any speedup?

    def get_seqs(self, indices, start=0, end=None):

        #Collect the sequences as a hash of arrays of strings.  Could use NumPy
        #for more efficient storage but it seems overkill as the number of reads
        #being extracted is comparatively small.
        #This also ensures that all the indices are ints
        seq_collector = { int(idx) : [] for idx in indices }

        #To build the sequence we have to loop over all the .bcl.gz files in the
        #cycle folders.  These are all named C[num].1 where num is 1-308 (unpadded).
        #The first base of the read will not be in C1.1 because there is an adapter,
        #and also for paired-end reads it doesn't make sense to read all bases as this
        #will run though onto the other end.  However, we can worry about this later.
        if end is None:
            end = self.num_cycles

        for cycle in range(start, end):
            cycle_dir = os.path.join(self.data_dir, 'C%i.1' % (cycle + 1))
            cycle_file = os.path.join(cycle_dir, '%s.bcl.gz' % self.prefix)

            fh = gzip.open(cycle_file, 'rb')
            bcl_header = fh.read(4)

            #The BCL header should be a fixed length depending on the machine type.
            #This assertion will fail when I try on output from different machines
            assert struct.unpack('<I', bcl_header)[0] == 4309650

            for idx in sorted(seq_collector.keys()):
                fh.seek(idx + 4)
                #Is reading bytes 1 at a time slow?  I'd imagine that internal cacheing
                #negates any need for chunked reads at this level.
                base_byte, = struct.unpack('B', fh.read(1))

                base = 'N'
                #qual = 0
                if base_byte:
                    base = ('A', 'C', 'G', 'T')[(base_byte & 0b11000000) >> 6]
                    #qual = (base_byte & 0b00111111) << 2 #Manual says "shifted by 2 bits" but I'm unsure

                seq_collector[idx].append(base)

        #TODO
        #Don't forget the accept/reject flag from the .filter file!

        #Remap the arrays into strings
        # return dict( idx : (nuc_string, flag) )
        return { idx : ( ''.join(seq), False ) for idx, seq in seq_collector.items() }

