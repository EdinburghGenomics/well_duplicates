#!/usr/bin/env python

import struct
import sys

fn = sys.argv[1]
sys.stderr.write(fn)
f = open(fn, 'rb')
o = open(sys.argv[2],'w')
sys.stderr.write(sys.argv[2])
try:
    buf = f.read(12)
    o.write("##%s\n"%(str(struct.unpack('=ifI', buf))))
    buf = f.read(8)
    while len(buf) == 8:
        t = struct.unpack('=ff', buf)
        x=int(round( 10 * t[0] + 1000))
        y=int(round( 10 * t[1] + 1000))
        o.write("%s\t%s\n"%(x,y))
        buf = f.read(8)
finally:
    o.close()
    f.close()
