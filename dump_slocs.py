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

f = None
try:
    f = open(sys.argv[1], 'rb')
except IndexError:
    f = sys.stdin

o = sys.stdout

try:
    buf = f.read(12)
    #First 3*4 bytes are the header
    header = struct.unpack('<ifI', buf)
    o.write("##%s\n"%(str(header)))
    buf = f.read(8)
    clusternum = 0
    clusternumpad = len(str(header[2]))
    while len(buf) == 8:
        # Each following 8 bytes are a co-ordinate pair as detailed in
        # https://broadinstitute.github.io/picard/javadoc/picard/picard/illumina/parser/readers/LocsFileReader.html
        # and
        # https://www.biostars.org/p/51681/
        t = struct.unpack('<ff', buf)
        x = int(t[0] * 10.0 + 1000.5)
        y = int(t[1] * 10.0 + 1000.5)
        #Write out "%i %i:%i\n"% (clusternum,x,y)
        #But I made it ugly by also interpolating the padding size
        o.write(("%%0%ii %%i:%%i\n"%clusternumpad) % (clusternum,x,y))
        clusternum += 1
        buf = f.read(8)
finally:
    o.close()
    f.close()
