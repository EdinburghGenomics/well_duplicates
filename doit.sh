#!/bin/bash -l
set -e ; set -u

### Note: This script is unlikely to be useful outside of Edinburgh Genomics ###

# This script aims to be something that can be executed as a cron job at, say,
# 15 minute intervals.  Therefore it needs to:
# 1) Run very fast if there is nothing to do - CHECK
# 2) Not get in a tizz if two instances start at once - CHECK
# 3) Not get stuck in a retry-loop if a job repeatedly fails - Handled in Snakefile
# 4) Push summary data to web1 where we can look upon it - CHECK
# 5) Add results to the Wiki as a sub-page of the run - Handled in Snakefile
# 6) Refuse to run on the backup headnode - CHECK
# 7) Log to a sensible location - CHECK

# Here's the fix for 6 and 2...
# Refuse to run on a machine other than headnode1
if [[ "${HOSTNAME%%.*}" != headnode1 ]] ; then
    echo "This script should only be run on headnode1"
    exit 1
fi

# 2) If $0 is not a canonical path, gripe
if [[ $(readlink -f "$0") != "$0" ]] ; then
    echo "You need to run this script by absolute path: $(readlink -f "$0")"
    exit 1
fi

# 3) Refuse to run twice (depends on 2 for reliable detection)
oldpid=`ps -lax | grep "^0.*/bash[ ]$0" | awk '{print $3}' | grep -vx $$` || true
if [[ -n "$oldpid" ]] ; then
    echo "$0 already running with PID $oldpid"
    exit 1
fi

# 4) This script already runs with "bash -l" in order to set up the SGE environment
#    and the extended PATH.  Otherwise I'd have to source /etc/profile.d/sge.sh and
#    add to the PATH here.

# Now we can send any further output to the log
LOGFILE=/ifs/runqc/WellDuplicates/logs/autorun.`date +%Y%m%d`.log
exec >>"$LOGFILE"

# Most things are dealt with by the Snakefile.  I need to run it over all the runs
# and push any new results that appear to web1.
# I also need to handle logging.
NOW=`date +%s`
SNAKEFILE="$(dirname $(readlink -f $0))"/Snakefile.count_and_push
export CLUSTER_QUEUE=casava

echo "=== Running at `date`. PID=$$, SNAKEFILE=$SNAKEFILE, CLUSTER_QUEUE=$CLUSTER_QUEUE ==="

for f in /ifs/seqdata/??????_[KE]00* ; do (
    echo "Trying to process $f"
    cd $f && "$SNAKEFILE" 2>&1 || true
) ; done

#This is not quoting-robust, but it will do
new_files=0
for f in `find /ifs/runqc/WellDuplicates/ -maxdepth 2 -mindepth 2 -name 10000targets_all_lanes.txt -newermt @$NOW` ; do
    new_files=1
    scp $f web1:/var/runinfo/WellDuplicates/$(basename $(dirname $f))_$(basename $f)
done

#And we can also make an overall summary table by munging the 10000targets_all_lanes.txt files.
if [ "$new_files" = 1 ] ; then
    ( echo Run $'\t'Lane{1..8}
      for f in /ifs/runqc/WellDuplicates/*/10000targets_all_lanes.txt ; do
          grep '^Level: 1' $f | sed 's/.*(\(.*\))$/\1/;H;$!d;g;s/\n/'$(basename $(dirname $f))'\t/;s/\n/\t/g'
      done
    ) > /ifs/runqc/WellDuplicates/summary.tsv
    set -x
    scp /ifs/runqc/WellDuplicates/summary.tsv web1:/var/runinfo/WellDuplicates/summary.tsv
fi

# That should doit.
echo "=== Finished run at `date`. PID=$$ ==="
