#!/urs/bin/env python3

import sys
import re
import io
import unittest
from unittest.mock import patch

try:
    from count_well_duplicates import output_writer, TALLY, LENGTH
except:
    #If this fails, you is probably running the tests wrongly
    print("****",
          "You want to run these tests from the top-level source folder by using:",
          "  python3 -m unittest test.test_count_well_duplicates",
          "or even",
          "  python3 -m unittest discover",
          "****",
          sep="\n")
    raise

#Ensure out assumptions about constants are right
assert TALLY == 0
assert LENGTH == 1

# dups_found = TILE_DUPL[tile][target][level][TALLY]
# targets_inspected = TILE_DUPL[tile][target][level][LENGTH]

# erm erm erm
# A lane has ~100 tiles.  A tile has maybe 2500 targets.
# A target has 5 levels (not including the centre).
# But here we have a lane with 1 tile, and the tile has 4 targets,
# and the target has 3 levels.
# See notebook sketch, which I will try to add add as a PNG.

# Each col = 1 level of target.  We're not recording lev 0!
TILE_DUPL = {'1208' : [
#           Level 1 T/L    2 T/L    3 T/L
                [ ( 0, 6), ( 0,12), ( 0,18) ],  #<-- 1 row = 1 target of this tile
                [ ( 2, 6), ( 1,10), ( 0,12) ],
                [ ( 3, 6), ( 1,10), ( 1,12) ],
                [ ( 0, 6), ( 1,12), ( 0,18) ] ] }

EXPECTED_OUT_1 = """
Lane:1   Tile: 1208   Targets: 4/4
Level:1   Wells: 24   Tot: 5   Hit: 2   AccO: 2   AccI: 3
Level:2   Wells: 44   Tot: 3   Hit: 3   AccO: 3   AccI: 3
Level:3   Wells: 60   Tot: 1   Hit: 1   AccO: 3   AccI: 1
LaneSummary:1   Targets: 4/4
Level:1   Wells:24   Tot:5 (0.208)   Hit:2 (0.500)   AccO:2 (0.500)   AccI:3 (0.750)
Level:2   Wells:44   Tot:3 (0.068)   Hit:3 (0.750)   AccO:3 (0.750)   AccI:3 (0.750)
Level:3   Wells:60   Tot:1 (0.017)   Hit:1 (0.250)   AccO:3 (0.750)   AccI:1 (0.250)
"""

BAD_TILE = {'1209':  [ ] }  #No valid targets.
BAD_TILE.update(TILE_DUPL)

# If we add the bad tile, what do we get?

EXPECTED_OUT_2 = """
Lane:1   Tile: 1208   Targets: 4/4
Level:1   Wells: 24   Tot: 5   Hit: 2   AccO: 2   AccI: 3
Level:2   Wells: 44   Tot: 3   Hit: 3   AccO: 3   AccI: 3
Level:3   Wells: 60   Tot: 1   Hit: 1   AccO: 3   AccI: 1
Lane:1   Tile: 1209   Targets: 0/4
Level:1   Wells: 0   Tot: 0   Hit: 0   AccO: 0   AccI: 0
Level:2   Wells: 0   Tot: 0   Hit: 0   AccO: 0   AccI: 0
Level:3   Wells: 0   Tot: 0   Hit: 0   AccO: 0   AccI: 0
LaneSummary:1   Targets: 4/8
Level:1   Wells:24   Tot:5 (0.208)   Hit:2 (0.500)   AccO:2 (0.500)   AccI:3 (0.750)
Level:2   Wells:44   Tot:3 (0.068)   Hit:3 (0.750)   AccO:3 (0.750)   AccI:3 (0.750)
Level:3   Wells:60   Tot:1 (0.017)   Hit:1 (0.250)   AccO:3 (0.750)   AccI:1 (0.250)
"""

class TestCountWellDuplicates(unittest.TestCase):

    #Capture sys.stdout - standard Mock procedure.
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_output_writer_1lane_full(self, mock_stdout):

        output_writer(1, TILE_DUPL)

        self._rescmp(mock_stdout, EXPECTED_OUT_1)

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_output_writer_badlane_full(self, mock_stdout):

        output_writer(1, BAD_TILE)

        self._rescmp(mock_stdout, EXPECTED_OUT_2)

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_output_writer_badlane_brief(self, mock_stdout):

        output_writer(1, BAD_TILE)

        self._rescmp(mock_stdout, EXPECTED_OUT_2, -4)

    def _rescmp(self, ioobj, astring, start=0, end=None):

        #The thing that got printed...
        lines1 = ioobj.getValue().split("\n")

        #The string we expected, ignoring leading newline
        #and swapping consecutive spaces for tabs.
        lines2 = [ re.sub('\s\s+', '\t', s) for s in astring.lstrip().split("\n") ]
        lines2 = lines2[start:end]

        #And now we can compare!
        self.assertEqual(lines1, lines2)
