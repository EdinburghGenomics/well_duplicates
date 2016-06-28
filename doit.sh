#!/bin/bash
set -e ; set -u

### Note: This script is unlikely to be useful outside of Edinburgh Genomics ###

# This script aims to be something that can be executed as a cron job at, say,
# 15 minute intervals.  Therefore it needs to:
# 1) Run very fast if there is nothing to do - CHECK
# 2) Not get in a tizz if two instances start at once - Handled in Snakefile
# 3) Not get stuck in a retry-loop if a job repeatedly fails - Handled in Snakefile
# 4) Push summary data to web1 where we can look upon it - CHECK
# 5) Add results to the Wiki as a sub-page of the run - Handled in Snakefile

# So most things are dealt with by the Snakefile.  I need to run it over all the runs
# and push any new results that appear.
# I also need to handle logging.
LOGFILE=/ifs/runqc/WellDuplicates/logs/autorun.`date +%Y%m%d`.log
exec >>"$LOGFILE"

NOW=`date +%s`
SNAKEFILE="$(dirname $(readlink -f $0))"/Snakefile.count_and_push
export CLUSTER_QUEUE=casava

echo "=== Running at `date`.  SNAKEFILE=$SNAKEFILE, CLUSTER_QUEUE=$CLUSTER_QUEUE ==="

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
    scp /ifs/runqc/WellDuplicates/summary.tsv web1:/var/runinfo/WellDuplicates/summary.tsv
fi

# That should doit.
