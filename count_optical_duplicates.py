#!/usr/bin/env python

from optparse import OptionParser
import sys
import logging
import os
from itertools import islice
import Levenshtein


def parse_coord_file(coord_file):
    try:
        coord_fh = open(coord_file, 'r')
    except IndexError:
        exit(1)

    coords = {}
    while True:
        coord_record = list(islice(coord_fh, 4))
        coords[coord_record[0]] = coord_record[1], coord_record[2], coord_record[3]
    coord_fh.close()

def sort_coords(coords):

    # unravel coords into one sorted list
    coord_list = []
    for key in coords:
        coord_list.append(key)
        coord_list.extend(coords[key][0])
        coord_list.extend(coords[key][1])
        coord_list.extend(coords[key][2])

    return coord_list.sort()

def get_edit_distance(str1, str2):

    return Levenshtein.distance(str1, str2)


def get_hamming_distance(str1, str2):

    return Levenshtein.hamming(str1, str2)


def main():

    sys.stderr.write("%s\n"%(str(sys.maxint)))
    # Setup options
    optparser = _prepare_optparser()
    (options, args) = optparser.parse_args()
    # verify options
    arg_pass = _verify_option(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)





def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-f coord_file> [-e edit_distance -n sample_size -l level]"""
    description = """This Script creates or executes commands that will assess optical duplicates
    in a run without mapping. Reads within level l of a selected reads from the coordinate file will be assessed for
    Levenshtein (edit) distance.
    """

    prog_version = "0.1"
    optparser = OptionParser(version=prog_version, description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h", "--help", action="help", help="show this help message and exit.")
    optparser.add_option("-f", "--coord_file", dest="coord_file", type="string",
                         help="The file containing the random sample per tile.")
    optparser.add_option("-e", "--edit_distance", dest="edit_distance", type="int",
                         help="max edit distance between two reads to count as duplicate")
    optparser.add_option("-n", "--sample_size", dest="sample_size", type="int", default=10000,
                         help="number of reads to be tested for optical duplicates (max number of prepared clusters is 10000 at the moment)")
    optparser.add_option("-l", "--level", dest="level", type="int", default=3,
                         help="levels around central spot to test, max = 3")

    return optparser


def _verify_option(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    if not options.coord_file:
        logging.error("You must specify a coordinates file.")
        arg_pass = False
    return arg_pass


if __name__ == "__main__":
    main()
