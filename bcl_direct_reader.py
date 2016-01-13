#!/usr/bin/python

from __future__ import print_function, division, absolute_import

"""
A module to grab sequence reads direct from the .bcl and .filter files
outputted by Illumina.  The motivation is that for some QC tasks we want
to gran a small subsample of reads, and getting these from the FASTQ is
very inefficient, especially once they are demultiplexed.

We assume not only the BCL format but also the standard directory
layout as specified in
https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf
"""

class BCLReader(object):

    def __init__(self, location="."):
        
    def get_seq_by_loc(lane, swath, tile, index, start=0, end=None):
        
        return nuc_string, flag

    def load_tile(lane, swath, tile):

        return tile # all tile data in memory

    #Option 3
    def get_tile(lane, swath, tile):

        return tile_ptr

class Tile:

    def get_seqs( [indices], start=0, end=None ):

        return dict( idx => (nuc_string, flag) )
