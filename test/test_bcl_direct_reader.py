#!python
from __future__ import print_function, division, absolute_import

import sys
import unittest

try:
    from bcl_direct_reader import BCLReader
except:
    #If this fails, you is probably running the tests wrongly
    print("****",
          "You want to run these tests by using:",
          "  python -m unittest test.test_bcl_direct_reader",
          "or even",
          "  python -m unittest discover",
          "****",
          sep="\n")
    raise

# Note that these tests are dependent on specific files being present in
# /ifs/seqdata/150715_K00169_0016_BH3FGFBBXX

TEST_PROJ = '/ifs/seqdata/150715_K00169_0016_BH3FGFBBXX'

class TestBCLReader(unittest.TestCase):

    def test_invalid_project(self):
        #What if we open the wrong folder?
        self.assertRaises(OSError, BCLReader, '/not/a/folder')

        #Or point to the Data subfolder in a valid project?
        self.assertRaises(OSError, BCLReader, TEST_PROJ + '/Data')

    def test_open_project(self):

        proj = BCLReader(TEST_PROJ)

        #For our test project we should see 8 lanes
        self.assertEqual(len(proj.lanes), 8)

    def test_get_seq_by_loc(self):
        
        #For this project, we can look at /ifs/seqarchive/150715_K00169_0016_BH3FGFBBXX/150715_K00169_0016_BH3FGFBBXX_1_1.sanfastq.gz
        #and see that the first sequence is:
        # @K00169:16:H3FGFBBXX:1:1101:1162:1791 1:N:0:ANNTCA
        # AATAGTCAGGTTAAATTTAATGTGACNNNTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGANTNAGNNNNN...

        # So the cluster co-ordinates are (1162:1791).  Running
        # $ dump_slocs.py s.locs | grep '1162:1791'
        # Gives us:
        # 0070657 1162:1791  (counting clusters from 0, not 1)

        proj = BCLReader(TEST_PROJ)

        seq, accept = proj.get_seq_by_loc(1, 1101, 70657, end=30)

        self.assertEqual(seq, 'AATAGTCAGGTTAAATTTAATGTGACNNNTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGANTNAGNNNNN')
        self.assertEqual(accept, true)
