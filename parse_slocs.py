#!/usr/bin/env python

import struct
import math


fn = 'X_s.locs'
f= open(fn, 'rb')
o = open('X_s.locs.txt','w')

try:
    buf = f.read(12)
    o.write("##%s\n"%(str(struct.unpack('=ifI', buf))))
    buf = f.read(8)
    while len(buf) == 8:
        t = struct.unpack('=ff', buf)
        x = int(t[0] * 10.0 + 1000.5)
        y = int(t[1] * 10.0 + 1000.5)
        yield (x,y)
        buf = f.read(8)
finally:
    o.close()
    f.close()

def yield_coords(o):