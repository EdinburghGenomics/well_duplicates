#!/usr/bin/env python

"""

input: sample_size n, sequencer_type (hiseq4000, hiseqx), levels (max 3)
return: dictionary of surrounding cluster indexes for n randomly selected wells
"""

import random
import sys
#from dump_slocs import yield_coords
import count_optical_duplicates
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


# maximum pixel distence between wells at a given level, required for edges of the flow cell
LEVEL_1_MAX_DIST = 22
LEVEL_2_MAX_DIST = 42
LEVEL_3_MAX_DIST = 62

DEF_SEED = 13
DEF_SEQ = "hiseq_4000"

DEF_SAMPLE_SIZE = 10000
MAX_CLUSTERS_4000 = 4312395
MAX_CLUSTERS_X = 6470949


def get_random_array(r_max, r_l, seed):
    if seed:
        random.seed(seed)
    ra = [random.randrange(0, r_max) for i in range(r_l)]
    return ra


def get_indexes(cluster_x, cluster_y, slocs_fh):

    l1_index = []
    l2_index = []
    l3_index = []

    # reset slocs file handle to position 12 (i.e. after the header)
    slocs_fh.seek(12)
    for coords in enumerate(yield_coords(slocs_fh)):
        (x,y) = coords[1]
        cluster_index = coords[0]
        dist = count_optical_duplicates.get_distance(cluster_x,cluster_y,x,y)
        if dist <=LEVEL_1_MAX_DIST:
            l1_index.append(cluster_index)
        elif dist <= LEVEL_2_MAX_DIST:
            l2_index.append(cluster_index)
        elif dist <= LEVEL_3_MAX_DIST:
            l3_index.append(cluster_index)

    return l1_index,l2_index,l3_index


def _prepare_argparser():
    """Prepare argparser object. New options will be added in this
    function first.
    """
    description = """This Script creates a list of n random cluster coordinates and
    the the indexes of surrounding coordinates with distances 1-3 from a HiSeq 4000 or HiSeqX s.locs file.
    """

    prog_version = "0.1"
    parser = ArgumentParser(description=description, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--slocs", dest="slocs", type=str,
                         help="The slocs file to analyse.")
    parser.add_argument("-t", "--type", dest="stype", type=str, default=DEF_SEQ,
                         help="the sequencer type, hiseq_4000 or hiseq_x")
    parser.add_argument("-s", "--seed", dest="seed", type=int, default=None,
                         help="Seed for the random read selection")
    parser.add_argument("-n", "--sample_size", dest="sample_size", type=int, default=DEF_SAMPLE_SIZE,
                         help="number of n random clusters")

    return parser.parse_args()

def _verify_option(args):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    if not args.slocs:
        sys.stderr.write("You must specify a read file.")
        arg_pass = False
    return arg_pass

def yield_coords(f):
    #You must have read the header first
    assert f.tell() >= 12

    buf = f.read(8)
    clusternum = 0
    while len(buf) == 8:
        # Each following 8 bytes are a co-ordinate pair as detailed in
        # https://broadinstitute.github.io/picard/javadoc/picard/picard/illumina/parser/readers/LocsFileReader.html
        # and
        # https://www.biostars.org/p/51681/
        t = struct.unpack('<ff', buf)
        x = int(t[0] * 10.0 + 1000.5)
        y = int(t[1] * 10.0 + 1000.5)

        yield (x,y)
        clusternum += 1
        buf = f.read(8)

def main():

    args = _prepare_argparser()

    # verify options
    arg_pass = _verify_option(args)
    if not arg_pass:

        sys.stderr.write("Non valid arguments: exit")
        sys.exit(1)

    # generate random list depending on MAX_CLUSTERS and sample_size
    # default sequencer is hiseq_4000
    random_sample = []
    if args.type == DEF_SEQ:
        random_sample = get_random_array(MAX_CLUSTERS_4000, args.sample_size, args.seed)
    elif args.type == "hiseq_x":
        random_sample = get_random_array(MAX_CLUSTERS_X, args.sample_size, args.seed)
    else:
        sys.stderr.write(args.type + ": Illigal argument, has to be hiseq_4000 or hiseq_x")

    slocs_fh = None
    try:
        slocs_fh = open(args.slocs, 'rb')
    except IndexError:
        slocs_fh = sys.stdin

    coord_dict = {}

    for coord in random_sample:
        slocs_fh.seek(12+(coord*8))  # 12 bytes for header, 8 byte per record, counting starts at 0
        cluster_x, cluster_y = yield_coords(slocs_fh)
        l1, l2, l3 = get_indexes(cluster_x, cluster_y, slocs_fh)
        coord_dict[coord] = l1, l2, l3

    for key in coord_dict.keys():
        print key + '\n' + '\n'.join(coord_dict[key]) + '\n'


main()