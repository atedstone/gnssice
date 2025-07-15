# gnss.sh

export GNSS_WORK='/scratch/hislide-processing/'

# Level 0 to Level 1 processing (only)
export GNSS_PATH_RINEX_DAILY=${GNSS_WORK}rinex_daily
export GNSS_PATH_RINEX_OVERLAP=${GNSS_WORK}rinex_overlap
export GNSS_PATH_RAWDATA=${GNSS_WORK}raw
export GNSS_PATH_TRACK_OUT=${GNSS_WORK}processed_track
export GNSS_RINEX_OBSERVER='SNF FlowState'
export GNSS_RINEX_INSTITUTION='University of Lausanne'

# Level 1 and Level 2 final output locations
export GNSS_L1DIR='/scratch/hislide-level1/'
export GNSS_L2DIR='/scratch/hislide-level2/'
