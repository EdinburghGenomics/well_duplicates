#!/usr/bin/env python

"""

input: sample_size n, sequencer_type (hiseq4000, hiseqx), levels (max 3)
return: dictionary of surrounding cluster indexes for n randomly selected wells
"""

import random
import sys
import struct
import math
# from dump_slocs import yield_coords
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# maximum pixel distence between wells at a given level, required for edges of the flow cell
LEVEL_1_MAX_DIST = 22
LEVEL_2_MAX_DIST = 42
LEVEL_3_MAX_DIST = 62

DEF_SEED = 13

DEF_SAMPLE_SIZE = 2500


def get_random_array(r_max, r_l, seed):
    if seed:
        random.seed(seed)
    ra = random.sample(range(r_max),r_l)
    return ra

def get_distance(x1, y1, x2, y2):

    dist = math.sqrt((x2-x1)**2+(y2-y1)**2)
    return dist


def get_indexes(cluster_coord, cluster_x, cluster_y, slocs_fh):
    """

    :rtype: object
    """
    l1_index = []
    l2_index = []
    l3_index = []

    # reset slocs file handle to position 12 bytes (i.e. after the header) or 5000*8bytes (0 or 5000 lines) before cluster_coord
    # TODO max distance for hiseq 4000, need to check for X
    offset = max([12, 12 + (cluster_coord - 5000) * 8])
    offset_coord = max([0, cluster_coord - 5000])
    sys.stderr.write("%s\n"%offset)
    slocs_fh.seek(offset)
    for coords in enumerate(yield_coords(slocs_fh)):
        (x, y) = coords[1]
        cluster_index = coords[0] + offset_coord
        dist = get_distance(cluster_x, cluster_y, x, y)
        if 1 < dist <= LEVEL_1_MAX_DIST:
            l1_index.append(cluster_index)
        elif LEVEL_1_MAX_DIST < dist <= LEVEL_2_MAX_DIST:
            l2_index.append(cluster_index)
        elif LEVEL_2_MAX_DIST < dist <= LEVEL_3_MAX_DIST:
            l3_index.append(cluster_index)
        if cluster_index > cluster_coord + 5000:
            break
    return l1_index, l2_index, l3_index


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
        sys.stderr.write("You must specify a read file.\n")
        arg_pass = False
    return arg_pass


def yield_coords(f):
    # You must have read the header first
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

        yield (x, y)
        clusternum += 1
        buf = f.read(8)


def main():
    args = _prepare_argparser()

    # verify options
    arg_pass = _verify_option(args)
    if not arg_pass:
        sys.stderr.write("Non valid arguments: exit")
        sys.exit(1)

    # TODO some logging, but needs a logger implementation
    sys.stderr.write("seed: %s\n" % (args.seed))
    sys.stderr.write("sample size: %s\n" % (args.sample_size))

    slocs_fh = None
    try:
        slocs_fh = open(args.slocs, 'rb')
    except IndexError:
        slocs_fh = sys.stdin

    # get MAX_CLUSTERS from header of s.locs file
    buf = slocs_fh.read(12)
    header = struct.unpack('=ifI', buf)
    MAX_CLUSTERS = int(header[2])
    sys.stderr.write("Maximum number of cluster according to s.locs: %s\n"%MAX_CLUSTERS)

    # generate random list depending on MAX_CLUSTERS and sample_size
    # default sequencer is hiseq_4000

    random_sample = get_random_array(MAX_CLUSTERS, args.sample_size, args.seed)

    sys.stderr.write("%s\n" % (random_sample))


    coord_dict = {}

    for coord in random_sample:
        slocs_fh.seek(12 + (coord * 8))  # 12 bytes for header, 8 byte per record, counting starts at 0
        # TODO should be imported from dump_slocs or suchlike, but I only need the one line here
        buf = slocs_fh.read(8)
        t = struct.unpack('=ff', buf)
        cluster_x = int(t[0] * 10.0 + 1000.5)
        cluster_y = int(t[1] * 10.0 + 1000.5)

        l1, l2, l3 = get_indexes(coord, cluster_x, cluster_y, slocs_fh)
        assert coord not in coord_dict
        coord_dict[coord] = l1, l2, l3

    for key in coord_dict.keys():
        print "%s\n%s\n%s\n%s" % (key,
                                  ",".join(str(n) for n in coord_dict[key][0]),
                                  ",".join(str(n) for n in coord_dict[key][1]),
                                  ",".join(str(n) for n in coord_dict[key][2]))

    slocs_fh.close()


main()
