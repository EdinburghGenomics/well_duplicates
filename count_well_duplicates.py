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
TALLY  = 0
LENGTH = 1

#Used for sequences returned from BCL Direct Reader
SEQUENCE  = 0
QUAL_FLAG = 1

def get_edit_distance(str1, str2):
    return Levenshtein.distance(str1, str2)


def get_hamming_distance(str1, str2):
    return Levenshtein.hamming(str1, str2)

def log(msg):
    print(str(msg), file=sys.stderr)

def output_writer(lane, sample_size, lane_dupl, levels=0, verbose=False):
    """ Reports on the lane by totting up the values in lane_dupl.
        The lane and sample_size arguments are added to the report but
        are not used in any calculations.

        If you want to understand what this is actually doing, look at
        the test code in test_count_well_duplicates.py

        Note it would be possible to output the data as it is generated
        without gethering the whole lot in RAM, but this seems a pointless
        optimization.
    """

    #First infer the number of levels, if not provided explicitly
    if not levels:
        for atile in lane_dupl.values():
            #Could be a duff tile?
            if len(atile) > 0:
                #No, it's OK.  Use the length of the first target.
                levels = len(atile[0])
                break

    #Grand totals...
    # Targets that got sampled (ie. had a valid read at the centre
    tot_targets = 0
    # Wells that got examined, at each level
    tot_wells = [0] * levels
    # Dups found at each level (total wells)
    tot_dups = [0] * levels
    # Dups found at each level (counting 1 per target per level)
    tot_hits = [0] * levels
    # Accumulated hits counting from the inside out
    tot_acco = [0] * levels
    # and counting from the outside in
    tot_acci = [0] * levels

    for tile in sorted(lane_dupl.keys()):

        tile_counts = lane_dupl[tile]
        targets = len(tile_counts)
        tot_targets += targets

        if verbose:
            print("Lane: %s\tTile: %s\tTargets: %i/%i" % (
                         lane,     tile,        targets,sample_size))

        #AccO and AccI require an explicit loop over targets.
        #I could tally the other things in this loop too but to me it makes
        #the code less readable.
        #Not that the code is pretty, but test coverage assures me it's good.
        acco = [0] * levels
        acci = [0] * levels
        for targ in tile_counts:
            seen_hit = 0
            for lev in range(levels):
                if targ[lev][TALLY]:
                    seen_hit = 1
                acco[lev] += seen_hit
            seen_hit = 0
            for lev in reversed(range(levels)):
                if targ[lev][TALLY]:
                    seen_hit = 1
                acci[lev] += seen_hit

        for lev in range(levels):

            wells = sum(targ[lev][LENGTH] for targ in tile_counts)
            dups = sum(targ[lev][TALLY] for targ in tile_counts)
            hits = sum(bool(targ[lev][TALLY]) for targ in tile_counts)

            if verbose:
                print("Level: %i\tWells: %i\tDups: %i\tHit: %i\tAccO: %i\tAccI: %i" % (
                              lev+1,     wells,    dups,    hits,     acco[lev],acci[lev]))

            tot_wells[lev] += wells
            tot_dups[lev] += dups
            tot_hits[lev] += hits

            tot_acco[lev] += acco[lev]
            tot_acci[lev] += acci[lev]

    #And report
    print("LaneSummary: %s\tTiles: %i\tTargets: %i/%i" % (
                        lane,      len(lane_dupl),
                                                tot_targets,
                                                   sample_size*len(lane_dupl) ))

    for lev in range(levels):
        print("Level: %i\tWells: %i\tDups: %i (%.3f)\t" % (
                      lev+1,     tot_wells[lev],
                                           tot_dups[lev],
                                               tot_dups[lev] / tot_wells[lev]) +
              "Hit: %i (%.3f)\tAccO: %i (%.3f)\tAccI: %i (%.3f)" % (
                    tot_hits[lev],
                        tot_hits[lev] / tot_targets,
                                     tot_acco[lev],
                                         tot_acco[lev] / tot_targets,
                                                      tot_acci[lev],
                                                          tot_acci[lev] / tot_targets)
             )

def main():
    # Setup options
    args = parse_args()

    if args.quiet:
        global log
        log = lambda *args: None

    lanes = args.lane.split(',') if args.lane else range(1, 8+1)

    max_tile = 24 #Works for Highseq X
    if args.stype == HIGHSEQ_4000:
        max_tile = 28

    tiles = []
    for swath in [11, 12, 21, 22]:
        for tile in range(1,max_tile+1):
            tiles.append("%s%02d" % (swath, tile))

    #If tiles are specified check that all are valid.
    if args.tile_id:
        for t in args.tile_id.split(','):
            assert t in tiles, "%s is not a valid tile for a %s" % (t, args.stype)
        tiles = args.tile_id.split(',')

    targets = load_targets( filename = args.coord_file,
                            levels = args.level+1,
                            limit = args.sample_size)
    bcl_reader = bcl_direct_reader.BCLReader(args.run)

    for lane in lanes:

        lane_dupl = {}
        for tile in tiles:
            log("Reading tile %s in lane %s" % (tile, lane))
            tile_bcl = bcl_reader.get_tile(lane, tile)

            #This actually reads the sequence data from the BCL into RAM
            seq_obj = tile_bcl.get_seqs(targets.get_all_indices(), args.start, args.end)
            log("Got %i sequences" % len(seq_obj))

            #Each entry in lane_dupl dict is a list of valid (ie. centre seq passed QC)
            #targets for this tile.
            lane_dupl[tile] = []

            for target in targets:

                center = target.get_centre()
                #log("Center: %s"%center)

                # if the center sequence does not pass the pass filter we don't assess edit distance
                # as large number of Ns compared to other reads with large number of Ns results in
                # small edit distance
                if not seq_obj[center][QUAL_FLAG]:
                    continue
                center_seq = seq_obj[center][SEQUENCE]

                #Add a placeholder for the new stats
                target_stats = [None] * args.level
                lane_dupl[tile].append(target_stats)

                for level in range(args.level):
                    #The level variable now runs from 0, but the target levels run from
                    #1 because 0 is the centre, so be careful!
                    dups = 0
                    well_indices = list(target.get_indices(level+1))
                    assert len(well_indices) > 0
                    for well_index in well_indices:
                        well_seq = seq_obj[well_index][SEQUENCE]
                        dist = get_edit_distance(center_seq, well_seq)

                        if dist <= args.edit_distance:
                            dups += 1
                            log("Center seq: %s" % center_seq)
                            log("well seq: %s" % well_seq)
                            log("edit distance: %s" % dist)

                    #Save a tuple of (TALLY, LENGTH)
                    target_stats[level] = (dups, len(well_indices))

            #log(lane_dupl)
        #Write output per lane
        output_writer(lane, len(targets), lane_dupl, verbose = not args.summary_only)


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
                        help="levels around central spot to test, max = 5")
    parser.add_argument("-s", "--stype", dest="stype", required=True, choices={HIGHSEQ_4000, HIGHSEQ_X},
                        help="Sequencer model")
    parser.add_argument("-r", "--run", dest="run", required=True,
                        help="path to base of run, i.e /ifs/seqdata/150715_K00169_0016_BH3FGFBBXX")
    parser.add_argument("-t", "--tile", dest="tile_id", type=str,
                        help="comma-separated list of specific tiles on a lane to analyse," +
                             " four digits, follow Illumina tile numbering")
    parser.add_argument("-i", "--lane", dest="lane", type=str,
                        help="comma-separated list of specific lanes to analyse, 1-8")
    parser.add_argument("-x", "--start", dest="start", type=int, default=50,
                        help="Starting base position for the slice of read to be examined")
    parser.add_argument("-y", "--end", dest="end", type=int, default=100,
                        help="Final base position for the slice of read to be examined")
    parser.add_argument("-S", "--summary-only", action="store_true",
                        help="Only print the summary per lane, not for every tile")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="No log output")
    parser.add_argument("--version", action="version", version=str(__VERSION__))

    return parser.parse_args()

if __name__ == "__main__":
    main()
