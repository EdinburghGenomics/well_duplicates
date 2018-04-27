#!/usr/bin/env python3
"""

input: sample_size n, sequencer_type (hiseq4000, hiseqx), levels (max 3)
return: dictionary of surrounding cluster indexes for n randomly selected wells
"""
__AUTHORS__ = ['Judith Risse', 'Tim Booth']
__VERSION__ = 0.2

import random
import sys
import struct
import math
# from dump_slocs import yield_coords
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# maximum pixel distence between wells at a given level, required for edges of the flow cell
# going out past 5 steps doesn't work properly!
MAX_DISTS = [1, 22, 42, 62, 82, 102]

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


def get_indexes(cluster_coord, cluster_x, cluster_y, slocs_fh, levels=5):
    """

    :rtype: object
    """
    MAX_SEARCH_AREA = 20000

    # I tried replacing this with "[[]] * levels" but that yields a list of N
    # references to a single list.  Oops.
    l_index = [[] for l in range(levels)]

    # reset slocs file handle to position 12 bytes (i.e. after the header)i
    # or 5000*8bytes (0 or 5000 lines) before cluster_coord
    # TODO max distance for hiseq 4000, need to check for X
    offset_coord = max([0, cluster_coord - MAX_SEARCH_AREA])
    offset = offset_coord * 8 + 12 #Byte offset in the file
    log(offset)
    slocs_fh.seek(offset)
    for coords in enumerate(yield_coords(slocs_fh)):
        (x, y) = coords[1]
        cluster_index = coords[0] + offset_coord
        dist = get_distance(cluster_x, cluster_y, x, y)

        for lev in range(levels): # 0,1,2,3,4
            if MAX_DISTS[lev] < dist <= MAX_DISTS[lev+1]:
                l_index[lev].append(cluster_index)

        #Stop searching when we get over 20000 records away
        if cluster_index > cluster_coord + MAX_SEARCH_AREA:
            break

    #Ensure we got something at every level
    for lev, wells in enumerate(l_index):
        if not l_index[lev]:
            raise RuntimeError(
                "Got no wells for cluster %s at (%s,%s) level %s",
                                         (cluster_coord,
                                                 cluster_x,
                                                    cluster_y,lev) )

    return l_index


def parse_args():
    """Prepare argparser object. New options will be added in this
    function first.
    """
    description = """This Script creates a list of n random cluster coordinates and
    the the indexes of surrounding coordinates with distances 1-3 from a HiSeq 4000 or HiSeqX s.locs file.
    """

    parser = ArgumentParser(description=description, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--slocs", dest="slocs", type=str, required=True,
                        help="The slocs file to analyse.")
    parser.add_argument("-s", "--seed", dest="seed", type=int, default=None,
                        help="Seed for the random read selection")
    parser.add_argument("-n", "--sample_size", dest="sample_size", type=int, default=DEF_SAMPLE_SIZE,
                        help="number of n random clusters")

    return parser.parse_args()

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

def log(msg):
    print(str(msg), file=sys.stderr)

def main():
    args = parse_args()

    # TODO some logging, but needs a logger implementation
    log("seed: %s" % (args.seed))
    log("sample size: %s" % (args.sample_size))

    slocs_fh = None
    try:
        slocs_fh = open(args.slocs, 'rb')
    except IndexError:
        slocs_fh = sys.stdin

    # get MAX_CLUSTERS from header of s.locs file
    buf = slocs_fh.read(12)
    header = struct.unpack('=ifI', buf)
    MAX_CLUSTERS = int(header[2])
    log("Maximum number of cluster according to s.locs: %s" % MAX_CLUSTERS)

    # generate random list depending on MAX_CLUSTERS and sample_size
    # default sequencer is hiseq_4000

    random_sample = get_random_array(MAX_CLUSTERS, args.sample_size, args.seed)

    log(random_sample)


    coord_dict = {}

    for coord in random_sample:
        slocs_fh.seek(12 + (coord * 8))  # 12 bytes for header, 8 byte per record, counting starts at 0
        # TODO should be imported from dump_slocs or suchlike, but I only need the one line here
        buf = slocs_fh.read(8)
        t = struct.unpack('=ff', buf)
        cluster_x = int(t[0] * 10.0 + 1000.5)
        cluster_y = int(t[1] * 10.0 + 1000.5)

        all_levs = get_indexes(coord, cluster_x, cluster_y, slocs_fh)
        assert coord not in coord_dict
        coord_dict[coord] = all_levs

    for key in coord_dict.keys():
        #Print the key on a line followed by a comma-separated list of
        #coords for each level out on one line each.
        print(str(key))
        for l in coord_dict[key]:
            print(",".join( map(str,l) ))

    slocs_fh.close()


main()
