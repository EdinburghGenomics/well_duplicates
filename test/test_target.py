#!python
from __future__ import print_function, division, absolute_import

import sys
import unittest
import time

try:
    # Adding this to sys.path helps the test work if you just run it directly.
    sys.path.insert(0,'.')
    from target import load_targets
except:
    #If this fails, you is probably running the tests wrongly
    print("****",
          "You want to run these tests from the top-level source folder by using:",
          "  python -m unittest test.test_target",
          "or even",
          "  python -m unittest discover",
          "****",
          sep="\n")
    raise

# These tests use a subset of the file hiseq_4000_10000clusters.list

TEST_FILE = 'test/small.list'
BAD_FILE_1 = 'test/bad1.list'
BAD_FILE_2 = 'test/bad2.list'

class TestTargetReader(unittest.TestCase):

    all_targets = None

    #Load the targets each time as the file is so small
    def setUp(self):
        self.all_targets = load_targets(TEST_FILE)

    def test_load_subset(self):
        sub_targets = load_targets(TEST_FILE, levels=2)

        self.assertEqual(sub_targets.levels, 2)

    def test_load_limit(self):
        lim_targets = load_targets(TEST_FILE, limit=2)

        self.assertEqual(len(lim_targets), 2)

        #Test iteration over all targets
        count = 0
        for targ in lim_targets:
            count +=1

        self.assertEqual(count, 2)

    def test_get_all_indices(self):
        #Use the first two targets, as above
        lim_targets = load_targets(TEST_FILE, levels=3, limit=2)

        self.assertEqual(
            set(lim_targets.get_all_indices(0)),
            set((1998850, 3178500))
        )

        self.assertEqual(
            set(lim_targets.get_all_indices(1)),
            set(map(int,(
                "1997278,1997279,1998849,1998851,2000420,2000421," +
                "3176929,3176930,3178499,3178501,3180071,3180072"
                ).split(",")))
        )

        self.assertEqual(
            set(lim_targets.get_all_indices(None)),
            set(map(int,(
                "1998850,1997278,1997279,1998849,1998851,2000420," +
                "2000421,1995707,1995708,1995709,1997277,1997280," +
                "1998848,1998852,2000419,2000422,2001991,2001992," +
                "2001993,3178500,3176929,3176930,3178499,3178501," +
                "3180071,3180072,3175357,3175358,3175359,3176928," +
                "3176931,3178498,3178502,3180070,3180073,3181641," +
                "3181642,3181643"
                ).split(",")))
        )

    def test_load_badfile(self):

        #bad1 has a blank line at the end
        self.assertRaises(ValueError, load_targets, BAD_FILE_1)

        #bad2 has the last line missing
        self.assertRaises(AssertionError, load_targets, BAD_FILE_2)

    def test_num_levels(self):
        all_targets = self.all_targets

        #I know there are 4 levels in the test file
        self.assertEqual(all_targets.levels, 4)

        self.assertEqual(all_targets.get_target_by_centre(196654).get_levels(), 4)

    def test_num_targets(self):
        all_targets = self.all_targets

        #There should be 7 of them
        self.assertEqual(len(all_targets), 7)

        #Ditto if we do it this way
        self.assertEqual(len(set(all_targets.get_all_indices(0))), 7)

    def test_lookups(self):
        all_targets = self.all_targets

        #If I look up point 196654 it should be the centre of cluster 6
        res = all_targets.get_from_index(196654)

        target_for_196654 = res[0][0]
        self.assertEqual(res, [ ( target_for_196654, 0 ) ])

        self.assertEqual(target_for_196654.get_centre(), 196654)

        res2 = target_for_196654.get_indices(1)
        self.assertEqual(res2, list(map(int,'195083,195084,196653,196655,198225,198226'.split(','))))

        #Getting all indices should be the same as getting them by level
        gathered_indices = set()
        for lev in range(4):
            gathered_indices.update(target_for_196654.get_indices(lev))
        self.assertEqual(gathered_indices, set(target_for_196654.get_indices()))

    def test_multiple_appearances(self):
        all_targets = self.all_targets

        #The file should contain exactly 213 distinct spots, as revealed by
        #sed 's/,/\n/g' test/small.list | sort -u | wc
        self.assertEqual(len(set(all_targets.get_all_indices())), 213)

        #Point 1030466 should pop up twice in rank 2 and once in rank 3
        res = all_targets.get_from_index(1030466)

        self.assertEqual(len(res), 3)

        self.assertEqual(sorted([ x[1] for x in res ]), [2,2,3])

    def test_bad_add(self):
        all_targets = self.all_targets

        self.assertRaises(Exception, all_targets.add_target, [(1, 2), (3, 4)])

        #This should complain about the number of levels
        self.assertRaises(AssertionError, all_targets.add_target, [(111,), (112, 113, 114, 115)])

        #This should work, but only once
        sub_targets = load_targets(TEST_FILE, 2)
        sub_targets.add_target( [(111,), (112, 113, 114, 115)] )
        self.assertRaises(Exception, sub_targets.add_target, [(111,), (112, 113, 114, 115)])

if __name__ == '__main__':
    unittest.main()
