#!python
from __future__ import print_function, division, absolute_import

import sys
import unittest
import time
import random

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

    def test_get_seq(self):
        """Tests fetching a single sequence
        """

        #For this project, we can look at /ifs/seqarchive/150715_K00169_0016_BH3FGFBBXX/150715_K00169_0016_BH3FGFBBXX_1_1.sanfastq.gz
        #and see that the first sequence is:
        # @K00169:16:H3FGFBBXX:1:1101:1162:1791 1:N:0:ANNTCA
        # AATAGTCAGGTTAAATTTAATGTGACNNNTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGANTNAGNNNNN...

        # So the cluster co-ordinates are (1162:1791).  Running
        # $ dump_slocs.py s.locs | grep '1162:1791'
        # Gives us:
        # 0070657 1162:1791  (counting clusters from 0, not 1)

        known_seq = 'AATAGTCAGGTTAAATTTAATGTGACNNNTTATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGANTNAGNNNNN'

        proj = BCLReader(TEST_PROJ)

        #For speed, just get the first 30 bases
        got_seq, got_accept = proj.get_seq(1, 1101, 70657, end=30)

        self.assertEqual(got_seq, known_seq[0:30])
        self.assertEqual(got_accept, True)

        #And test the start/end by looking at bases 20-40
        got_seq, got_accept = proj.get_seq(1, 1101, 70657, start=20, end=40)

        self.assertEqual(got_seq, known_seq[20:40])
        self.assertEqual(got_accept, True)

    #This takes about a minute the first time, and about 30 seconds on subsequent
    #reads from the same tile.  Comment out the next line to try it.
    @unittest.skip("slow test")
    def test_worst_case_speed(self):

        #The slowest single read is the full sequence of the very last cluster
        #on a tile, which on our test data is 4309650-1

        proj = BCLReader(TEST_PROJ)

        start_time = time.time()

        tile = proj.get_tile(1, 1101)
        slow_seq, _ = tile.get_seqs([4309649])[4309649]

        self.assertEqual(len(slow_seq), tile.num_cycles)

        print("\n*** Slow read took %.2f seconds." % (time.time() - start_time))

    @unittest.skip("slow test")
    def test_multi_read_speed(self):

        # Due to the way I read the data, getting bases from many spots should
        # be fairly efficient.
        proj = BCLReader(TEST_PROJ)
        tile = proj.get_tile(1, 1101)

        # How to generate some consistent lists?
        fetchlists = {}
        sample_sizes = [1,10,100,1000,10000,100000]

        start_time = time.time()
        for count in sample_sizes:
            random.seed(0)
            fetchlists[count] = random.sample(range(4309650), count)
        print("\n*** Making the number lists took %.2f seconds." % (time.time() - start_time))

        #Now see about fetching
        for count in sample_sizes:
            start_time = time.time()
            res = tile.get_seqs(fetchlists[count], start=20, end=40)
            print("\n*** Fetching 20 bases from %i sequences took %.2f seconds." %
                                               (count,           (time.time() - start_time))
                 )

        #In my mind, fetching the first 10000 seqs should be faster than fetching just the last
        #in the file, but maybe not...
        start_time = time.time()
        res = tile.get_seqs([4309649], start=20, end=40)
        print("\n*** Fetching 20 bases from sequence 4309649 took %.2f seconds." % (time.time() - start_time))

        start_time = time.time()
        res = tile.get_seqs(range(100000), start=20, end=40)
        print("\n*** Fetching 20 bases from the first 100000 seqs took %.2f seconds." % (time.time() - start_time))

    def test_invalid_get_seqs(self):

        #What if I ask for a sequence that's out of range?
        proj = BCLReader(TEST_PROJ)

        self.assertRaises(IndexError, proj.get_seq, 1, 1101, 5000000)

    def test_multiple_get_seqs(self):
        """Tests batch sequence fetching, and also tests the accept/reject flag
        """

        #The library was designed to extract multiple sequences at once, so test this.
        #Also test that the accept/reject flag is being got correctly.  Again I'll use
        #sequences from the file above, tile 1101
        good_seqs = ( ( 70657 , 'AATAGTCAGG' ), # at 1162:1791
                      ( 70658 , 'TGTGGCATTT' ), # at 1182:1791
                      ( 70660 , 'TCAGAATCAG' ), # at 1223:1791
                      ( 89638 , 'TTGCTTATCA' ), # at 4024:2002 (line 50001 in the .sanfasq)
                      ( 89639 , 'GCCTTATGGC' ), # at 4045:2002
                      ( 89641 , 'TACTGAGAAG' )) # at 4085:2002

        #I can't find these in the FASTQ so assumed they were rejected reads
        #The sequences were obtained from running bcl_direct_reader.py so unlike
        #the good_seqs above are only valid for regression testing.
        bad_seqs = ( ( 70659, 'TTGGACGAGG' ),
                     ( 89640, 'CCCNNNNNNN' ),
                     (     0, 'NNNNNNNNNN' ),
                     (     1, 'NNNNNNNNNN' ),
                     (     2, 'NNNNNNNNNN' ))

        #Synthesize the expected result:
        expected_result = { clus[0] : ( clus[1], True )
                            for clus in good_seqs  }

        expected_result.update( { clus[0] : (clus[1], False)
                                  for clus in bad_seqs } )

        #print(expected_result)

        proj = BCLReader(TEST_PROJ)
        tile = proj.get_tile(1, 1101)

        result = tile.get_seqs([clus[0] for clus in good_seqs + bad_seqs], end=10)

        self.assertEqual(result, expected_result)

        #TODO - if I implement any caching then test that I get the same result on a second call

