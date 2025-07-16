# Example/template for GNSSICE environment variables
# Suggest copying to home directory and then modifying to your system.

export GNSS_WORK='/work/atedstone/gnss_example/'

# Level 0 to Level 1 processing (only)
export GNSS_PATH_RAWDATA=${GNSS_WORK}raw
export GNSS_PATH_RINEX_DAILY=${GNSS_WORK}rinex_daily
export GNSS_PATH_RINEX_OVERLAP=${GNSS_WORK}rinex_overlap
export GNSS_PATH_IONEX_DAILY=${GNSS_WORK}ionex_daily
export GNSS_PATH_IONEX_OVERLAP=${GNSS_WORK}ionex_overlap
export GNSS_PATH_SP3_DAILY=${GNSS_WORK}sp3_daily
export GNSS_PATH_SP3_OVERLAP=${GNSS_WORK}sp3_overlap
export GNSS_PATH_TRACK_OUT=${GNSS_WORK}processed_track

export GNSS_RINEX_OBSERVER='SNF FlowState'
export GNSS_RINEX_INSTITUTION='University of Lausanne'


# Level 1 and Level 2 final output locations
export GNSS_L1DIR='/work/atedstone/flowstate-gnss-level1/'
export GNSS_L2DIR='/work/atedstone/flowstate-gnss-level2/'
