#!/bin/bash
# Tidy-up, produce files for archival for a given rover site.
# Don't need to treat base data as we no longer run our own. 

set -euo pipefail
# Turn on debugging
# set -x

# Note, this doesn't work with raw Leica files, but does work with raw GVTs
if [ $# -ne 2 ]; then
    echo "Usage: $0 <site> <year>" >&2
    exit 1
fi

site="$1"
year="$2"

# Generate Compact RINEX files
for f in rinex_daily/$site/*o; do
    gzf=${f:0:-1}d.gz
    if [ ! -f $gzf ]; then 
        # .yyo to .yyd
        rnx2crx $f
        # .gz
        gzip ${f:0:-1}d;
    else
        echo "Skipping $f (d.gz exists)"
    fi
done

# Zip the processing logs
zip track_processing_logs_${site}_${year}.zip processed_track/$site/*${year}_*.out
# And move them to a 'central' folder for easy transfer by scp
mkdir -p processing_logs
mv track_processing_logs_${site}_${year}.zip processing_logs/