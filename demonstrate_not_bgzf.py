#!python

# Run this with:
# env PYTHONPATH=$PYTHONPATH:/ifs/software/linux_x86_64/biopython/biopython-1.66/build/lib.linux-x86_64-2.7 \
#   /ifs/software/linux_x86_64/python/python_2.7.10_venv/bin/python demonstrate_not_bgzf.py

from __future__ import print_function, division, absolute_import

import sys, os
from Bio import bgzf

# The first should work.  The second should fail.
files_to_inspect = ['/ifs/software/linux_x86_64/biopython/biopython-1.66/Tests/GenBank/NC_000932.gb.bgz',
                    '/ifs/seqdata/150715_K00169_0016_BH3FGFBBXX/Data/Intensities/BaseCalls/L002/C1.1/s_2_1103.bcl.gz',
                   ]

#Plus anything on the command line you want to look at.
files_to_inspect.extend(sys.argv[1:])

for f in files_to_inspect:

    print("Inspecting %s" % os.path.basename(f))

    handle = open(f, "rb")

    try:
        for values in bgzf.BgzfBlocks(handle):
            print("  Raw start %i, raw length %i; data start %i, data length %i" % values)
    except Exception as e:
        print("  Error: " + str(e))

    handle.close()
