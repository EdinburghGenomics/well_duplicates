#!/bin/bash
set -euo pipefail
set -o noclobber

# This logic was broken out of Snakefile.count. It generates a targets file
# for a given s.locs file with a given number of targets. If the cache is available
# and an appropriate file already exists then the output will be a symlink.
# If the cache is not in use the output will be a regular file.
# A very basic locking mechanism ensures we don't write junk into the cache (a
# second process trying to write to the cache will just fail).

# Output will never be clobbered so if there's an old file you need to remove it
# first.
trap "echo 'Usage: get_cached_targets.sh <locs_file> <target_count> <output_file>'" EXIT
locs_file="$1"
target_count="$2"
output_file="$3"
trap - EXIT

# First ensure I'll run the corresponding well dups scripts
WD_ROOT="$(dirname $(readlink -f $0))"
PATH="$WD_ROOT:$PATH"

# Then work out where the cache directory is at
CLUSTER_LISTS="${CLUSTER_LISTS:-$WD_ROOT/cluster_lists}"

if [ -d "$CLUSTER_LISTS" ] ; then
    md5=`md5sum "$locs_file" | awk '{print $1}'`

    cached_list="$CLUSTER_LISTS/${target_count}clusters_${md5}.list"

    if ! [ -e "${cached_list}.done" ] ; then
        trap "rm '$cached_list'" EXIT
        prepare_cluster_indexes.py -n "$target_count" -f "$locs_file" > "$cached_list"
        touch "${cached_list}.done"
        trap - EXIT
    fi
    # Shall the link be relative?
    if [ "${cached_list:0:1}" = / ] ; then
        ln -s "$cached_list" "$output_file"
    else
        ln -sr "$cached_list" "$output_file"
    fi
else
    # No-cache mode it is, then
    prepare_cluster_indexes.py -n "$target_count" -f "$locs_file" > "$output_file"
fi
