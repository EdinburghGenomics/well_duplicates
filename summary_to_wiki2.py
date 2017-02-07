#!/usr/bin/env python3
### NOTE : This script is unlikely to be useful outside of Edinburgh Genomics. ###

import os, sys, re

"""Like summary_to_wiki.py.  Reads the 10000targets_all_lanes.txt summary text file
   as generated by the Snakefile and spits out just the fraction of repeats for
   each lane, in a 2-column table.
   Note the file is actually created by running 'tail' on the individual results
   which is where the ==> Headings <== comes from.
   I could make the code output Wiki markup directly but I don't want to have
   that in the scripts, so this is done in the style of old-school "perl -p"
   regex-powered munging.
   The output going to be added as comments so I need to generate HTML markup.
"""

lane = "0"

def munge(line):
    global lane

    sum_mo = re.match(r"LaneSummary: (\d+).*Tiles:", line)
    if sum_mo:
        lane = sum_mo.group(1)

    """
    lev1_mo = re.match("Level: 1\s", line)
    if lev1_mo:
        #headings = [i.split(": ")[0] for i in line.split("\t")]
        vals = [i.split(": ",1)[1] for i in line.split("\t")]

        frac_dup = re.search(r"\((.*)\)", vals[-1]).group(1)
        perc_dup = round(float(frac_dup) * 100, len(frac_dup))

        return trow(lane, "%s %%" % perc_dup)
    """
    #dup_mo = re.match(r"Overall duplication .*: ([0-9.%]+)", line)
    dup_mo = re.match(r"Picard-equivalent duplication score v2: *([0-9.%]+)", line)
    if dup_mo:
        #Put a space before the % sign
        perc_dup = dup_mo.group(1).replace("%", " %")
        return trow(lane, perc_dup)

#HTML silliness
def trow(*args):
    return "<tr>" + ''.join("<td>%s</td>" % a for a in args) + "</tr>"

def h3(*args):
    return ''.join("<h3>%s</h3>" % a for a in args)

def tstart(*headers):
    return '<table>\n<tr>' + ''.join("<th>%s</th>" % h for h in headers) + "</tr>"

def tend():
    return '</table>'

if __name__ == '__main__':
    #pipeline mode
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL)

    print(h3("Well Duplicates Summary"))
    print(tstart("Lane","Est. Duplication"))

    for line in (x.rstrip("\n") for x in sys.stdin):
        munged = munge(line)
        if munged is not None:
            print(munged)

    print(tend())
