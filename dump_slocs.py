#!/usr/bin/env python

from __future__ import print_function, division, absolute_import

import struct
import math
import sys

#Python normally complains about being killed by SIG_PIPE, but we just want
#exit gracefully if that happens (eg. if output is piped to head)
#(http://docs.python.org/library/signal.html)
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

def main():
    f = None
    try:
        f = open(sys.argv[1], 'rb')
    except IndexError:
        f = sys.stdin

    o = sys.stdout

    try:
        header = yield_header(f)
        o.write("##%s\n"%(str(header)))
        clusternumpad = len(str(header[2]))

        for coords in enumerate(yield_coords(f)):
            #Write out "%i %i:%i\n"% (clusternum,x,y)
            #But I made it ugly by also interpolating the padding size
            o.write(("%%0%ii %%i:%%i\n"%clusternumpad) % (coords[:1] + coords[1][:]))
            #o.write(str(coords) + "\n")

    finally:
        o.close()
        f.close()

def yield_header(f):
    #First 3*4 bytes are the header
    #Only the third element is really useful as this is the number of locations
    #that should be in the file.
    buf = f.read(12)
    return struct.unpack('<ifI', buf)

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

main()
