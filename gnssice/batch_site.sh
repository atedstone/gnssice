#!/bin/bash
# Create RINEX files and run double-differencing for a given site.
# Before usage for a given year, make sure to update the hard-coded
# specifications of different base station periods.
# A.T., 31.08.2026

set -euo pipefail
# Turn on debugging
set -x

if [ $# -ne 4 ]; then
    echo "Usage: $0 <site> <type>(G or R) <year> <start_doy>" >&2
    exit 1
fi

site="$1"
type="$2"
year="$3"
start_doy="$4"


# ------------------------------------------------------------------------------
## Create site RINEX files

# not for L1200 sites!!
#process_rinex $site $type -overlap


# ------------------------------------------------------------------------------
## A-priori coordinates from the specified first DOY to process

file="rinex/${site}/${site}${start_doy}0.${year:2:2}o"
 
if [ ! -f "$file" ]; then
    echo "Error: file '$file' not found" >&2
    exit 1
fi
 
# Find the line containing the label, then pull out the first three
# whitespace-delimited fields (the X, Y, Z values).
line=$(grep 'APPROX POSITION XYZ' "$file" | head -n 1)
 
if [ -z "$line" ]; then
    echo "Error: no line containing 'APPROX POSITION XYZ' found in '$file'" >&2
    exit 1
fi
 
x=$(echo "$line" | awk '{print $1}')
y=$(echo "$line" | awk '{print $2}')
z=$(echo "$line" | awk '{print $3}')
 
# Remove any hyphen/minus signs
x=${x//-/}
y=${y//-/}
z=${z//-/}
 
result="${x} ${y} ${z}"
 
echo "$result"


# ------------------------------------------------------------------------------
## Launch TRACK
# Modify this section for each year of processing according to base
# station availability

process_dgps klsq $site $year $start_doy 135 -ap $result --unsup

process_dgps lrhp $site $year 136 300 --unsup

process_dgps klsq $site $year 301 365 --unsup

process_dgps klsq $site ${year+1} 1 125 --unsup


# ------------------------------------------------------------------------------
## Concatenate batches
# Modify this section for each year of processing according to base
# station availability

concat_daily_geod klsq $site $year $start_doy 135

concat_daily_geod lrhp $site $year 136 300

concat_daily_geod klsq $site $year 301 365

concat_daily_geod klsq $site ${year+1} 1 125


 
