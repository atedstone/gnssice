#!/bin/bash
# Create RINEX files and run double-differencing for a given site.
# Before usage for a given year, make sure to update the hard-coded
# specifications of different base station periods.
# A.T., 31.08.2026

set -euo pipefail
# Turn on debugging
# set -x

if [ $# -ne 4 ]; then
    echo "Usage: $0 <site> <type>(G or R) <year> <start_doy>" >&2
    exit 1
fi

site="$1"
type="$2"
year="$3"
start_doy="$4"

## =============================================================================

# Reusable function to locate a file named ssssddd0.YYo given a site,
# year, and starting day-of-year. If the file for start_doy doesn't
# exist, it searches forward for the next available day.
#
# Usage:
#   find_obs_file <site> <year> <start_doy> [search_dir]
#
# On success: prints the found filepath to stdout and returns 0.
# On failure: prints an error to stderr and returns 1.
# Claude Sonnet 5 Medium, Andrew Tedstone, 
# "Extract and format coordinates values from text file", 01.09.2026
 
find_obs_file() {
    local site="$1"
    local year="$2"
    local start_doy="$3"
    local search_dir="${4:-.}"   # optional 4th arg; defaults to current dir
 
    local yy="${year: -2}"
    local doy ddd candidate
 
    for (( doy=10#$start_doy; doy<=366; doy++ )); do
        ddd=$(printf "%03d" "$doy")
        candidate="${search_dir}/${site}${ddd}0.${yy}o"
 
        if [ -f "$candidate" ]; then
            echo "$doy"
            return 0
        fi
    done
 
    echo "Error: no file found for site '$site' on or after DOY $start_doy in year $year (dir: $search_dir)" >&2
    return 1
}



# ------------------------------------------------------------------------------
## Create site RINEX files

# not for L1200 sites!!
process_rinex $site $type -overlap


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
## Launch TRACK and concatenate to batches
# Modify this section for each year of processing according to base
# station availability
# Use find_obs_file() to avoid missing RINEX files crashing the workflow.

doy=$(find_obs_file $site $year $start_doy "rinex/${site}/")
#process_dgps klsq $site $year $doy 135 -ap $result --unsup
conc_daily_geod klsq $site $year $doy 135

doy=$(find_obs_file $site $year 136 "rinex/${site}/")
process_dgps lrhp $site $year $doy 300 --unsup
conc_daily_geod lrhp $site $year $doy 300

doy=$(find_obs_file $site $year 301 "rinex/${site}/")
process_dgps klsq $site $year $doy 365 --unsup
conc_daily_geod klsq $site $year $doy 365

doy=$(find_obs_file $site  $((year+1)) 1 "rinex/${site}/")
process_dgps klsq $site $((year+1)) $doy 125 --unsup
conc_daily_geod klsq $site $((year+1)) $doy 125

