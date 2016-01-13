#!/usr/bin/env python

from optparse import OptionParser
import sys
import logging
import os
from itertools import islice
import random
import math
import Levenshtein
from qc_utils import utils_logging

SEED = 13
SAMPLE_SIZE = 10000


def get_file_size(file_name):
    return os.stat(file_name).st_size


def get_bytes_fq_record(fq_file_name):
    fq = open(fq_file_name, 'r')
    rec = ''.join(islice(fq, 4))
    byte_record = len(rec)
    fq.close()
    return byte_record


def get_records_fq_file(fq_file_name):

    return get_file_size(fq_file_name)/get_bytes_fq_record(fq_file_name)


def get_random_array(r_max, r_l, seed):
    if seed:
        random.seed(seed)
    ra = [random.randrange(1, r_max) for i in range(r_l)]
    return ra


def get_distance(x1, y1, x2, y2):

    dist = math.sqrt((x2-x1)**2+(y2-y1)**2)
    return dist


def parse_fq_header(header):
    p = header.split(' ')[0].split(':')
    return p[4], p[5], p[6]


def get_edit_distance(str1, str2):

    return Levenshtein.distance(str1, str2)


def get_hamming_distance(str1, str2):

    return Levenshtein.hamming(str1, str2)


def main():

    sys.stderr.write("%s\n"%(str(sys.maxint)))
    # initialize the logging
    utils_logging.init_logging(logging.INFO)
    # Setup options
    optparser = _prepare_optparser()
    (options, args) = optparser.parse_args()
    # verify options
    arg_pass = _verify_option(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    #if options.debug:
    #    utils_logging.init_logging(logging.DEBUG)

    fq_rec = get_records_fq_file(options.read_1_file)
    b_fq = get_bytes_fq_record(options.read_1_file)
    sys.stderr.write("records in file: %s\n"%(fq_rec))
    b_rec = get_bytes_fq_record(options.read_1_file)
    sys.stderr.write("bytes per record: %s\n"%(b_rec))
    rand_array = get_random_array(fq_rec, options.sample_size, options.seed)

    fq1 = open(options.read_1_file, 'r')
    #fq2 = None
    #if options.read_2_file:
    #    fq2 = open(options.read_2_file, 'r')

    for r in rand_array:
        sys.stderr.write("record: %s\n"%(r))
        #b_rec = numpy.int64(b_rec)
        byts = b_rec*(r-1)
        sys.stderr.write("bytes: %s\n"%(byts))
        sys.stderr.write("should be: 8143916874\n")
        fq1.seek(byts)
        r1_header1, r1_seq, r1_header2, r1_qual = islice(fq1, 4)
        sys.stderr.write("fastq1:\n%s%s%s%s"%(r1_header1, r1_seq, r1_header2, r1_qual))
        t1, x1, y1 = parse_fq_header(r1_header1)

        sys.stderr.write("UP\n")
        for i in range(1,50):
            fq1.seek(b_rec*(r-i))
            r2_header1, r2_seq, r2_header2, r2_qual = islice(fq1, 4)
            #sys.stderr.write("fastq2:\n%s%s%s%s"%(r2_header1, r2_seq, r2_header2, r2_qual))
            t2, x2, y2 = parse_fq_header(r2_header1)
            dist = 'NA'
            if t1 == t2:
                dist = get_distance(int(x1), int(y1), int(x2), int(y2))

            h_dist = get_hamming_distance(r1_seq, r2_seq)
            l_dist = get_edit_distance(r1_seq, r2_seq)

            #sys.stderr.write("Read 1: tile: %s\tx: %s\ty: %s\n"%(t1, x1, y1))
            #sys.stderr.write("Read 2: tile: %s\tx: %s\ty: %s\n"%(t2, x2, y2))
            sys.stderr.write("Pixel distance: %s\tHamming distance: %s\tLevenshtein distance: %s\n"%(dist, h_dist, l_dist))
            #sys.stderr.write("%s%s\n"%(r1_seq, r2_seq))

        sys.stderr.write("DOWN\n")
        for j in range(1,50):
            fq1.seek(b_rec*(r+j))
            r2_header1, r2_seq, r2_header2, r2_qual = islice(fq1, 4)
            #sys.stderr.write("fastq2:\n%s%s%s%s"%(r2_header1, r2_seq, r2_header2, r2_qual))
            t2, x2, y2 = parse_fq_header(r2_header1)
            dist = 'NA'
            if t1 == t2:
                dist = get_distance(int(x1), int(y1), int(x2), int(y2))

            h_dist = get_hamming_distance(r1_seq, r2_seq)
            l_dist = get_edit_distance(r1_seq, r2_seq)

            #sys.stderr.write("Read 1: tile: %s\tx: %s\ty: %s\n"%(t1, x1, y1))
            #sys.stderr.write("Read 2: tile: %s\tx: %s\ty: %s\n"%(t2, x2, y2))
            sys.stderr.write("Pixel distance: %s\tHamming distance: %s\tLevenshtein distance: %s\n"%(dist, h_dist, l_dist))
            #sys.stderr.write("%s%s\n"%(r1_seq, r2_seq))



def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog <-f read_1_file> [-r read_2_file -s seed -n sample_size -d distance]"""
    description = """This Script create or execute commands that will assess optical duplicates
    in a run without mapping. It will assess n reads per tile from read_1_file and read_2_file,
    using s as random seed. Reads within d pixel of a selected read will be assessed for
    Levenshtein (edit) distance.
    """

    prog_version = "0.1"
    optparser = OptionParser(version=prog_version, description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h", "--help", action="help", help="show this help message and exit.")
    optparser.add_option("-f", "--read_1_file", dest="read_1_file", type="string",
                         help="The read 1 file to analyse.")
    optparser.add_option("-r", "--read_2_file", dest="read_2_file", type="string",
                         help="The read 2 file to analyse.")
    optparser.add_option("-s", "--seed", dest="seed", type="int", default=None,
                         help="Seed for the random read selection. Default: %default")
    optparser.add_option("-n", "--sample_size", dest="sample_size", type="int", default=10000,
                         help="number of reads to be tested for optical duplicates")
    optparser.add_option("-d", "--distance", dest="distance", type="int", default=60,
                         help="Pixel distance to select neighbouring reads. Default: %default")

    return optparser


def _verify_option(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    if not options.read_1_file:
        logging.error("You must specify a read file.")
        arg_pass = False
    return arg_pass


if __name__ == "__main__":
    main()
