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

# Here's the quick fix for 6...
# Refuse to run on headnode2
if [[ "${HOSTNAME%%.*}" == headnode2 ]] ; then
    echo "This script should not be run on headnode2" >&2
    exit 1
fi

# Settings specific to old/new clusters. Note that this stand-alone script should soon
# be superceded on the new cluster by a QC stage incorporated into Illuminatus.
if [ -e /lustre/software ] ; then
    WORKDIR_ROOT="$HOME/WellDuplicates"
    SEQDATA=/lustre/seqdata
else
    WORKDIR_ROOT=/ifs/runqc/WellDuplicates
    SEQDATA=/ifs/seqdata
fi
LOGFILE="$WORKDIR_ROOT/logs/autorun.`date +%Y%m%d`.log"

# Previously the script checked 'ps' for other instances, but I'm switching to the
# recommended 'flock' mechanism.
FLOCK_FILE="${TMPDIR:-/tmp}/flock_$(readlink -f "$0" | md5sum | awk '{print $1}')"
if [ "${FLOCK_ON:-0}" = 0 ] ; then
    # echo "Locking exclusively on $FLOCK_FILE, PID=$$"
    (   flock -n 9 || exit 33
        export FLOCK_ON=9
        source "$0" "$@"
    ) 9>"$FLOCK_FILE" ; rc="$?"
    if [ "$rc" != 0 ] ; then
        if [ "$rc" = 33 ] ; then
            echo "Failed to gain shared lock, PID=$$" >> "$LOGFILE"
        else
            #This should trigger an e-mail to the cron manager
            echo "Script exited with error $rc" >&2
        fi
        exit "$rc"
    fi
    #Spawned copy ran, nothing more to do.
    #echo "Exiting unlocked script, PID=$$"
    exit 0
fi
#echo "Locked on $FLOCK_FILE, PID=$$"

# 4) This script already runs with "bash -l" in order to set up the SGE environment
#    and the extended PATH.  Otherwise I'd have to source /etc/profile.d/sge.sh and
#    add to the PATH here.

# Now we can send any further output to the log
exec >>"$LOGFILE"

# Most things are dealt with by the Snakefile.  I need to run it over all the runs
# and push any new results that appear to web1.
# I also need to handle logging.
SNAKEFILE="$(dirname $(readlink -f $0))"/Snakefile.count_and_push
export CLUSTER_QUEUE=casava

echo "=== Running at `date`. PID=$$, SNAKEFILE=$SNAKEFILE, CLUSTER_QUEUE=$CLUSTER_QUEUE ==="

for f in "$SEQDATA"/??????_[KE]00* ; do (
    echo "Trying to process $f"
    export WORKDIR="$WORKDIR_ROOT/`basename $f`"

    cd $f && "$SNAKEFILE" 2>&1 || true
) ; done

# Copying to web1 has been removed. See GIT on 21/2/17 for the old version.

# That should doit.
echo "=== Finished run at `date`. PID=$$ ==="
