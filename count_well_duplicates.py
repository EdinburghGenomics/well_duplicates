#!/usr/bin/env python

from __future__ import division, print_function, absolute_import

__AUTHORS__ = ['Judith Risse', 'Tim Booth']
__VERSION__ = 0.2

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
from itertools import islice
import Levenshtein
import bcl_direct_reader
from target import load_targets

DEF_SEQ = HIGHSEQ_4000 = "hiseq_4000"
HIGHSEQ_X = "hiseq_x"

#Used for indexing tuple per target:
TALLY = 0
LENGTH = 1

def get_edit_distance(str1, str2):
    return Levenshtein.distance(str1, str2)


def get_hamming_distance(str1, str2):
    return Levenshtein.hamming(str1, str2)

def log(msg):
    print(str(msg), file=sys.stderr)

def output_writer(lane, tile_dupl, levels=0):
    To be written, in accordance with the test!

def main():
    # Setup options
    args = parse_args()

    lanes = args.lane.split(',') if args.lane else range(1, 8+1)

    tiles = []
    max_tile = 0
    if args.stype == HIGHSEQ_4000:
        max_tile = 28
    else:
        max_tile = 24

    if args.tile_id:
        tiles = args.tile_id.split(',')
    else:
        for swath in [11, 12, 21, 22]:
            for tile in range(1,max_tile+1):  # should be 24 for hiseq_X
                tile_id = "%s%02d" % (swath, tile)
                tiles.append(tile_id)

    targets = load_targets(args.coord_file, args.level+1)
    bcl_reader = bcl_direct_reader.BCLReader(args.run)

    for lane in lanes:
        # dict[tile]
        tile_dupl = {}
        for tile in tiles:
            tile_bcl = bcl_reader.get_tile(lane, tile)
            seq_obj = tile_bcl.get_seqs(targets.get_all_indices(), args.start, args.end)

            # Initialise tally and length counters for this tile
            tile_dupl[tile] = [{'tally': 0, 'length': 0} for level in range(0, args.level+1)]

            for target in islice(targets.get_all_targets(), 0, args.sample_size):

                center = target.get_centre()
                #log("Center: %s"%center)

                # if the center sequence does not pass the pass filter we don't assess edit distance
                # as large number of Ns compared to other reads with large number of Ns results in
                # small edit distance
                if not seq_obj[center][1]:
                    continue
                center_seq = seq_obj[center][0]

                for level in range(1, args.level+1):
                    l_dupl = []
                    assert(target.get_levels()>= level)
                    for well_index in target.get_indices(level):
                        well_seq = seq_obj[well_index][0]
                        dist = get_edit_distance(center_seq, well_seq)

                        if dist <= args.edit_distance:
                            l_dupl.append(1)
                            log("Center seq: %s" % center_seq)
                            log("well seq: %s" % well_seq)
                            log("edit distance: %s" % dist)
                        else:
                            l_dupl.append(0)
                    if sum(l_dupl) >= 1:
                        tile_dupl[tile][level]['tally'] += 1
                    tile_dupl[tile][level]['length'] += 1

                    # TODO - see if the lines above were really what you wanted and if
                    # so fix the comment in output_writer to reflect what we are actually
                    # recording.

            log(tile_dupl)
        output_writer(lane, tile_dupl, args.level)


def parse_args():
    description = """This script creates or executes commands that will assess well duplicates
    in a run without mapping. Reads within level l of a selected reads from the coordinate file
    will be assessed for Levenshtein (edit) distance.
    """

    parser = ArgumentParser(description=description, formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f", "--coord_file", dest="coord_file", required=True,
                        help="The file containing the random sample per tile.")
    parser.add_argument("-e", "--edit_distance", dest="edit_distance", type=int, default=2,
                        help="max edit distance between two reads to count as duplicate")
    parser.add_argument("-n", "--sample_size", dest="sample_size", type=int, default=2500,
                        help="number of reads to be tested for well duplicates (max number" +
                             " of prepared clusters is 10000 at the moment)")
    parser.add_argument("-l", "--level", dest="level", type=int, default=3,
                        help="levels around central spot to test, max = 3")
    parser.add_argument("-s", "--stype", dest="stype", required=True, choices={HIGHSEQ_4000, HIGHSEQ_X},
                        help="Sequencer model")
    parser.add_argument("-r", "--run", dest="run", required=True,
                        help="path to base of run, i.e /ifs/seqdata/150715_K00169_0016_BH3FGFBBXX")
    parser.add_argument("-t", "--tile", dest="tile_id", type=str,
                        help="specific tile on a lane to analyse, four digits, follow Illumina tile numbering")
    parser.add_argument("-i", "--lane", dest="lane", type=str,
                        help="specific lane to analyse, 1-8")
    parser.add_argument("-x", "--start", dest="start", type=int, default=50,
                        help="Starting base position for the slice of read to be examined")
    parser.add_argument("-y", "--end", dest="end", type=int, default=100,
                        help="Final base position for the slice of read to be examined")
    parser.add_argument("--version", action="version", version=str(__VERSION__))

    return parser.parse_args()

if __name__ == "__main__":
    main()
