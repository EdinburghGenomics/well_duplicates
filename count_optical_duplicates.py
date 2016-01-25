#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
import logging
import Levenshtein
import bcl_direct_reader
from target import load_targets


def get_edit_distance(str1, str2):
    return Levenshtein.distance(str1, str2)


def get_hamming_distance(str1, str2):
    return Levenshtein.hamming(str1, str2)


def main():
    # Setup options
    optparser = _prepare_argparser()
    args = optparser.parse_args()
    # verify options
    arg_pass = _verify_option(args)
    if not arg_pass:
        logging.critical("Non valid arguments: exit")
        sys.exit(1)

    lanes = range(1,8)
    if args.lane:
        lanes = [args.lane]

    tiles = []
    if args.tile:
        tiles = [args.tile_id]
    else:
        for swath in [11,12,21,22]:
            for tile in range(1,28):
                tile_id = "%s%s"%(swath, tile)
                tiles.append(tile_id)

    targets = load_targets(args.coord_file)
    bcl_reader = bcl_direct_reader.BCLReader(args.run)

    for lane in lanes:
        for tile in tiles:

            tile_bcl = bcl_reader.get_tile(lane, tile)

            seq_obj = tile_bcl.get_seqs(targets.get_all_indices())

            for target in targets.get_all_targets():
                center = target.get_centre()
                center_seq = seq_obj[center]
                for level in range(1,args.level):
                    l_dupl = []
                    for well_index in target.get_indices(1):
                        well_seq = seq_obj[well_index]
                        dist = get_edit_distance(center_seq,well_seq)
                        if dist <= args.edit_distance:
                            l_dupl.append(1)
                        else:
                            l_dupl.append(0)


def _prepare_argparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-f coord_file> [-e edit_distance -n sample_size -l level]"""
    description = """This Script creates or executes commands that will assess optical duplicates
    in a run without mapping. Reads within level l of a selected reads from the coordinate file will be assessed for
    Levenshtein (edit) distance.
    """

    prog_version = "0.1"
    parser = ArgumentParser(description=description, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--coord_file", dest="coord_file", type=str,
                         help="The file containing the random sample per tile.")
    parser.add_argument("-e", "--edit_distance", dest="edit_distance", type=int,
                         help="max edit distance between two reads to count as duplicate")
    parser.add_argument("-n", "--sample_size", dest="sample_size", type=int, default=2500,
                         help="number of reads to be tested for optical duplicates (max number of prepared clusters is 10000 at the moment)")
    parser.add_argument("-l", "--level", dest="level", type=int, default=3,
                         help="levels around central spot to test, max = 3")
    parser.add_argument("-r", "--run", dest="run", type=str,
                         help="path to base of run, i.e /ifs/seqdata/150715_K00169_0016_BH3FGFBBXX")
    parser.add_argument("-t", "--tile", dest="tile_id", type=str,
                         help="specific tile on a lane to analyse, four digits, follow Illumina tile numbering")
    parser.add_argument("-i", "--lane", dest="lane", type=int,
                         help="specific lane to analyse, 1-8")

    return parser


def _verify_option(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    if not options.coord_file:
        logging.error("You must specify a coordinates file.")
        arg_pass = False
    if not options.run:
        logging.error("You must specify a run folder")
        arg_pass = False
    return arg_pass


if __name__ == "__main__":
    main()
