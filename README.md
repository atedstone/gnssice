# GNSSIce - Processing of differential GNSS observations from ice sheets and glaciers

This package provides functionality to use TRACK to kinematically process differential GNSS/GPS observations collected on moving ice, and to post-process them to yield displacements and velocities.

If this is your first time using this package to process GNSS data, suggest reading the HISLIDE data user guide before continuing further.


## Sticky - Level-0 to Level-1 processing only: octopus.unil.ch quick-start

Setup two terminals. Use the first terminal for processing. 

Terminal 1:

```bash
# Connect to login node
ssh octo
# Connect to compute node
ssh [node]
# Set up environment
source gnss_process.sh
```

Where `gnss_process.sh` contains lines similar to this:

```bash
# Python env
mamba activate gnss
# TRACK
module load gamit_globk
# Change to working location
cd /work/atedstone/
# List contents
ls
# Set up lightweight server so we can see PNGs
python -m http.server 8080 &
```

In a second terminal, set up port forwarding to the lightweight HTTP server so that we can look at PNG files of outputs when needed.

Terminal 2:

```bash
ssh octo -N -L localhost:8080:[node]:8080
firefox localhost:8080
```

Once you are finished with processing, in the first terminal remember to stop the web server process:

```bash
# Find the process ID
ps -u
# Stop it
kill -9 <PID>
```


N.b. this quick-start assumes the following pre-requisites:

* Local machine is Unix or Mac, octopus is configured in SSH config with X11 forwarding enabled
* gnssice installed on octopus



## Overview of workflow

### Data levels

- Level-0 : Refers to both (i) raw binary data straight from the receivers (e.g. `*.ubx`) and (ii) the representation of those same data as RINEX files.
- Level-1 : GEOD data which have been kinematically post-processed with TRACK.
- Level-2 : Derived values, e.g. along- and across-track displacements.

### Processing Level-0 to Level-1

If you only want to process level-1 to level-2 then you can ignore this section.

1. Ensure everything is installed.
2. Set up the working environment with VMF, IONEX, SP3 files.
4. Copy base and rover files to working directory.
5. (Obtain third party data to cover gaps in our base record.)
6. `process_rinex`: Convert raw files to daily RINEX files (base and rover), and overlap them in windows if needed (usually best practise to do the overlapping.)
7. `process_dgps`: Do TRACK kinematic processing.
    6a. Process 'temporary' fixes taken during redrilling/flights.
8. `conc_daily_geod`: Concatenate daily TRACK GEOD files to year.


Example using site **lev5** with a Leica 1200+ receiver, 2021, days 129 to 242:

```bash
# First edit gnss.sh to set your desired paths, then source it
source gnss.sh 
# Prepare working directories, download VMF3, IONEX, SP3, overlap them.
setup_wd 2021 129 242 -overlap

# Copy the relevant TRACK cmd file template to your working directory
cp -r /location/<TRACK_cmd_filename> $GNSS_WORK
# Verify the contents of the cmd file, making sure that dcb.dat, antenna files are available.
vi <TRACK_cmd_filename>

# Make the raw data available (repeat for the base)
cd $GNSS_PATH_RAWDATA
mkdir lev5; cd lev5
cp -r /location/myrawdata .

# Process Leica raw files into RINEX (also do this for the base)
process_rinex lev5 LS -filename Default...m00 s -start 2021-05-09 -finish 2021-08-30 -overlap

# (If this receiver were a Cryologger GVT, we'd do this instead:)
process_rinex lev5 G -overlap

# Run TRACK
# First get a-priori coordinates, and increase sigmas for cmd file site sigmas
process_dgps rusb lev5 2021 129 129 -ap x y z
# Reduce site sigmas in cmd file then continue
process_dgps rusb lev5 2021 130 242 

# Post-process the batch of data you just processed.
# This creates parquet files in the processed_track/lev5/ directory.
# They have filenames rover_base_year_startDOY_endDOY_GEOD.parquet
conc_daily_geod rusb lev5 2021 129 242

# Eventually, combine batches of data together 
# This produces filename rover_startyear_startDOY_endyear_endDOY_geod.parquet
export_level1 lev5
```


### Processing Level-1 to Level-2

9. `calculate_local_origin`: Only if this is a new site, to calculate local origin position. If this is an existing site then you should already have a file 'origin_<site>.csv' available.
10. `gnss_disp_vel.py`: Preferentially run interactively/as notebook. Transform coordinates, filter data, convert to along/across-TRACK displacements, calculate velocities. Exclude periods of data based on the user-input exclusion file. N.b. use of this script requires care depending on the length of baseline and the speed of the site, check the script for more details!
11. Be sure to retain the post-processing ancillary files if they are to be used to process another batch of data from a site in the future. (rotation matrix, site origin)
12. `seasonal_annual_disp.py` : To calculate seasonal and annual displacements.

Continuing our example above: 

```bash
# Run next line only if no origin.csv file:
calculate_local_origin lev5 lev5_rusb_2021_129_2021_242_geod.parquet
# Calculate velocities of this batch of data...
gnss_disp_vel.py lev5 lev5_rusb_2021_129_2021_242_geod.parquet
# ...Or 'add' them onto an existing dataset.
gnss_disp_vel.py lev5  -f ...
```

## Installation

### The `gnssice` package

Set up a conda or mamba environment:

```bash
conda create -n gnss pip
conda activate gnss
```

If you haven't already got it, clone the gnssice repo, change into its directory then install in-place using pip (this will also install all dependencies):

```bash
pip install -e .
```

This means that it won't matter where you then do your main working, the scripts will always be accessible. *(Note: the command line interface uses entry points: see the full list in `setup.cfg`.)*

**Important**: Installing `gnssice` allows you to process Level-1 to Level-2 only. For Level-0 to Level-1, several other dependencies are required -- see below for information.


### Dependencies needed for Level-0 to Level-1 data

#### TRACK

This package wraps TRACK, part of GAMIT/GLOBK: http://geoweb.mit.edu/gg/. Usage requires a license, which needs to be requested from MIT in the case that your Institution is not already a User. Users are provided with an access password for GLOBK/Gamit/TRACK downloads.

You should check that the Gamit tables are sufficiently up to date for your captured GPS epochs. If they are not, update the GLOBK/Gamit installation following their instructions.


#### RINEX utilities

For the following, download the binaries relevant for your system and place them in your `bin`:

- `gfzrnx`, the GFZ Potsdam RINEX toolbox: https://gnss.git-pages.gfz-potsdam.de/gfzrnx/welcome/. 
- RINEX hatanaka compression tools: http://terras.gsi.go.jp/ja/crx2rnx.html
- For older receivers (e.g. Leica 1200), we still rely on TEQC, even though this is discontinued since 2019. Download an executable suitable for your system from http://facility.unavco.org/software/teqc/.
- For newer receivers (e.g. ublox), RTKLIB (https://www.rtklib.com/) is often suitable.  Download the latest source (https://github.com/tomojitakasu/RTKLIB) and then build the `convbin` utility (this package does not use any other RTKLIB utilities). For build instructions see Section 4 of the RTKLIB manual.


## Detailed Usage: Level-0 to Level-1

### Multi-year field campaigns

**This workflow treats data in 'batches'.** A batch of data corresponds to a specific base-rover combination over a specific time frame. Batches of data cannot span multiple years. So, if picking data from the field only once a year e.g. in springtime, the processing still needs to be split into two batches separated by the change in year.

**Therefore, set up a work environment dedicated to a single calendar year.**

### Summary of file types

In the following,

* `p`: provider (e.g. igs)
* `y`: year (`yy` corresponds to the two-digit year, e.g. 21 for 2021)
* `d`: day of year

* Orbit files: `pppyyddd0.sp3` where `p` is provider (e.g. igs), `y` is year, `d` is day of year
* RINEX files:
	- **!!! Note that we ONLY use the RINEX2 filename format, even for RINEX3 files !!!**
	- If compressed: `*.yyd`
	- If uncompressed: `*.yyo`
* TRACK command files: `TRACK_<base>_<rcvtype>.cmd`
* Files generated by TRACK:
	- Log file: `<rover>_<base>_<year>_<doy>.out`
	- GEOD results file: `<rover>_<base>_<year>_<doy>_GEOD.dat` 
	- NEU results file: `<rover>_<base>_<year>_<doy>_GEOD.dat` 
    - Figures: `<rover>_<base>_<year>_<doy>_NEU.png` 
	- Processing session log: `gps.TRACK.<rover>.log`


### Explanation of working directory and structure

This package relies on environment variables to find input and output file locations. You should copy and adapt `gnss.sh` to suit your needs.

The root of all processing is defined by `$GNSS_WORK`. All Level-0 to Level-1 processing commands will operate here, even if your terminal prompt location is elsewhere on the server.

Looking at the rest of `gnss.sh`, see that we distinguish between 'daily' versus 'overlapped' input files, each of which are stored in their own distinct folders. The 'overlap' folders will only exist if you explicitly request overlapping files.

Whether you use overlapping files or not, TRACK always expects the following locations inside the `$GNSS_WORK` directory to be available: 

- `vmf` 
- `sp3` (symlink)
- `ionex` (symlink)
- `raw/<site>/`
- `rinex/<site>/` (symlink)
- `processed_track/<site>/` 

If you use `setup_wd` (next section) then all the folders in the `gnss.sh` file and in the list above will get created automatically, apart from `processed_track` which gets made later by `process_dgps`. As indicated in the list above, some of these are actually symlinks to the underlying 'daily' or 'overlap' folders, depending on whether you are using overlapping files in your workflow.


### Preparations

**Easy preparations:** `source gnss.sh` then `setup_wd` to setup your working directory according to your environment variables, including creating all the folders needed to get started with analysis. It will download all SP3, IONEX and VMF3 files for your period of interest, and overlap them if requested. (This could take several tens of minutes if the period of interest is a whole year). If you subsequently need to extend the time series of these files, just run the command again with the updated start and end days.

**To overlap or not?** Short answer = yes, overlap. There is usually no good reason not to. This workflow will handle making overlapping files for you if you ask it to (RINEX, SP3 and IONEX files).

Following is an explanation of how do these preparations manually.

#### Get the orbit files

```bash
get_orbits <year> <start doy> <end doy>
```

Add `-overlap` if running the workflow with overlapping daily RINEX files.
	 
N.b. don't attempt this over a change in year, e.g. start of 360 and end of 4. Instead do two calls. The only problem will then be that days 365 and 1 don't contain the next and previous days data respectively. Just copy and paste from the relevant files into the next, or on unix command line use cat, e.g:

```bash
cat igs364.sp3 igs365.sp3 igs001.sp3 > igs365.sp3
```
   

#### Optional: Get 3rd party base RINEX files.

If required, obtain RINEX files from another site to cover the gaps.

Sometimes 3rd party base files can be downloaded directly from public archives using the GAMIT/GLOBK script `sh_get_rinex`: 

```bash
sh_get_rinex -archive sopac -yr 2011 -doy 0 -ndays 250 -sites kely
```

If files have `*.<yy>o` suffix you're all set, otherwise, if they are zipped, unzip the compressed rinex files using 7zip or whatever. Then convert to normal rinex i.e. from `*.10d` to `*.10o`, using CRX2RNX.

```bash	
for f in *; do
crx2rnx $f;
done
```	    
 
Overlap/window RINEX files using the `R` file type argument to `process_rinex`:

```bash
process_rinex site R -overlap
```


#### IONEX files

Recommended for long baselines, e.g. > 100 km. See also the TRACK help info and Doyle et al. (2014) supplementary methods.

There is a simple helper script, `get_ionex`, which can be modified to download IONEX files as needed.  

TRACK needs to be told about these files by adding the `IONEX_FILE` option in the cmd file:

```bash
 IONEX_FILE ionex/igsg<day>0.<S09>i
``` 


#### TRACK COMMAND file

Ensure that the TRACK command file is correct for the site you are processing. Basically, each combination of base station and receiver type needs its own command file. There are several details to get right.

* The filename must be in the format `track_<base>_<rcvtype>.cmd`, where:
	* `<base>` is the 4-letter base station identifier, e.g. rusb, klsq
	* `<rcvtype>` refers to the relevant entry in the receivers lookup table `TRACK_CMD_RCVR` in `gnss.py`. E.g. for a cryologger GVT use the value `gvt`.
* Check that the `site_pos` for the base station is correct.
* Ensure that the stated DCB file is available.
* Ensure that the `ante_off` settings are correct 
	- Check that the antenna set in the command file matches the antenna indicated in a RINEX file header for the site.
	- Check that the receiver type (last argument of `ante_off`) is correct
* Check that the stated antmod_file is available.
		
	
#### Convert the raw files to RINEX

Copy the raw files for the rover and base into the scratch space.

We can convert to RINEX files which are windowed to 28hrs duration: from 22:00 the previous day to 02:00 the next day. 

In your terminal window, make sure you're in your `$GNSS_WORK` directory.

Use `process_rinex`, run with `-h` to find out the options.

Make sure you create RINEX files for each site (i.e. run the script for each site).

After, you can delete any empty files:

	find -size 0c -delete


### Do the kinematic processing

Pre-requisites:

- A number of additional arguments are hard-set in the cmd files created for each base station, (e.g. levb, kely), e.g. site_stats and bf_set. These probably don't need to be modified from their default values...
- Open `process_dgps.py` and ensure that the default parameters for processing with your base station are set up - follow the format in the file. Again, set initial values based on recommendation above, unless you've more specific information to go on... 


#### A-priori coordinates

For first day for each site use a priori coordinates of each site derived from online 

- Often the best option: webapp.geod.nrcan.gc.ca/geod/tools-outils/ppp.php
- Maybe: http://apps.gdgps.net/apps_file_upload.php

Upload a rinex file and specify the day from which to return a priori coordinates.
Subsequently select no. The program takes the APR coordinates from the previous days results.

*In 2012 I had problems with coordinates from the above website. However, teqc estimates the daily positions and puts them in the top of daily rinex files - these seem to do the trick, so I used these instead. You might also find that the site_stats (see later) initially have to be loosened to deal with the APR coordinates. (AJT, September 2012).*


#### Processing

With loose `site_stats` in the TRACK cmd file, run only the first day of the site using a-priori coordinates:

```bash
process_dgps <base> <rover> <year> <start DOY> <start DOY> -ap <X> <Y> <Z>
```

Once this has completed, tighten the `site_stats` in the TRACK cmd file, then run:

```bash
process_dgps <base> <rover> <year> <start DOY> <end DOY>
```
	
You don't have to process an entire site in one session. Enter the start day and guesstimate an appropriate end day. If you get fed up before the end day is reached, do `CTRL+C` to break out/halt the process (preferably between processing two days of data, rather than during).

If providing a-priori coordinates (-ap): TRACK will not transfer negative values to the cmd file. Instead, put them positive here and make sure there is a negative sign (-) in the relevant field of the cmd file. Keep this negative sign in place for all subsequent processing of a site with a negative coordinate.

`process_dgps` has defaults for `ion_stats`, `MW_WL` and `LG` set. These vary depending on the specified base station: at the time of writing, levb and kely are both supported.

By default, TRACK is set up to accept the day's results automatically, using an RMS approach. (A Spearman correlation approach is also available). The thresholds are set in the arguments to gnss.Kinematic.TRACK().

TRACK can also be set up to require the user to accept every day manually - provide `use_auto_qa=False` as an argument to the call to `gnss.Kinematic.TRACK` in `process_dgps`. 

Initially, TRACK will try to use the default ion_stats, MW_WL and LG parameters as specified in `process_dgps.py`.

If you choose to reject TRACK's initial results based on the default parameters, you can enter your own, explanations below are based on the TRACK help file:

- `jump`, the largest jump in the ion delay allowed before bias flag introduced. This parameter sets `@ion_stats <jump>` in the TRACK cmd file. (String number: `<S06>`)
- `WL_Fact`, the MW-WL weighting. Weight to be given to deviation of MW-WL from zero. Default is 1 (ie., equal weight with LC residuals). Setting the value smaller will downweight the contribution of the MW-WL. For noisy or systematic range data (can be tested with a P1,P2 or PC solution), the WL_fact may be reduced.
	+ MW WL stands for Melbourne-Wubbena (MW) wide lane (this is the combination of range and phase that gives the difference between the L1 and L2 biases and it independent of the ionospheric delay and the geometry - TRACK help file line 266). The parameter is taken by the float_type command, setting `@FLOAT_TYPE <WL_Fact>` (string number: `<S04>`).
- `Ion_fact`, the LG combination weighting. Weight to be given to deviation of the Ionospheric delay from zero.  Default is 1 (i.e., ionospheric delay is assumed to be zero and given unit weighting inz deterimining how well a set of integer ambiquities fit the data.  On long baselines, value should be reduced.  Sets `<Ion_fact>` argument for @FLOAT_TYPE (string number `<S05>`).
	-  For long baselines (>20 km) this value should be reduced to give less weight to the ionospheric delay constraint.  For 100 km baselines, 0.1 seems to work well.  With very good range data (ie., WL ambiquities all near integer values), this factor can be reduced.

Particularly if processing long baselines, processing times of >15 mins per day are not unusual.		

If results fail the automated quality check, they will pop up in a figure window, and the RMS values will appear in the terminal window. The figure window will block input to the terminal window until it is closed.

You'll be asked to input some information for logging - a quality identifier, and an optional longer comment. These are appended the the log file.

Quality identifiers:

* G: Good
* O: Ok
* B: Bad
* A: Accepted Automatically (by RMS threshold or Spearman test)

If RMS high and results not good, try (in the following order):

1. Increase `ion_stats`, up to maybe 3 or 4. Choose the combination that gives the lowest rms and produces the closest to a straight line when plotted.
2. Reduce `WL_fact` from default=1 to say 0.5 if data appear noisy.
3. Possibly try reducing `Ion_fact` - especially for sites further away where the ionosphere has an increasingly large effect.

In addition, sometimes increasing `site_stats` in the TRACK cmd file can help when results are of poor quality. But, the defaults set in the cmd file reflect our predictions that sites are unlikely to be > 10 metres away from a priori position (ish), and should not have moved more than 2 metres compared to the last epoch. So, increasing these values shouldn't really help... (but it does sometimes!)

Possible errors:

- `SP3 Interpolation Errors`. This may mean that there are problems with one or more satellites in the SP3 orbit files. Open the .out file for the day which just failed and go to the bottom of it. Check the line(s) which say "Interpolation error PRN" - the number directly after this is the satellite. Open the cmd file for the relevant base station and add: `exclude_svs <satellite_number>`
- (e.g.) `TRACK IOSTAT error:  IOSTAT error      2 occurred opening gps_config/TRACK_kely.cmd in MAKEXK` - this probably means you're trying to process while not in your working directory.
- `!!APR .LC file does not exist! Terminating processing.` If you are expecting a .LC file to exist this may again mean that you are trying to process while not in your working directory.
- `STOP FATAL Error: Stop from report_stat.` Have you remembered to update the RINEX file suffix in the TRACK CMD files to reflect the year you are processing?
- Something about fortran files: you probably haven't loaded gamit module - see start of this doc.
- `For        0 Double differences: Average RMS        NaN mm` and lots of clock warnings in the TRACK .out file: check that your a-priori coordinates for the rover site are good. If you took them directly from the header of the first RINEX file then there a higher likelihood that they are bad; either look at other RINEX files or consider getting a PPP solution.

Old ideas on trying to improve results - sometimes work but we're even less certain why!

- Increase site_stats in the .cmd (TRACK_kely.cmd/TRACK_levb.cmd) file (df=10 10 10 1 1 1, try 100 100 100 10 10 10 first, then gradually increase)
- Increase bf_set in the .cmd file (bf= 2 40, try 5 100 and then increase gradually)
- Increase or decrease ion_stats in the input parameters (df=1, try 0.3 first then 2,3 if that dosnt work)
	
	
If a day or more accidentally or otherwise get missed during processing and you need to go back to re-process them:

- Find the NEU and GEOD files for the day prior that on which you need to re-start processing.
- Make sure they are in your workspace directory.
- Rename to `TRACK.NEU.<site>.LC` and `TRACK.GEOD.<site>.LC`
- This ensures that TRACK uses the previous day's coordinates for a priori estimate positioning.
- If you haven't finished processing the rest of the season, be sure to change filenames back to the appropriate LC files.


#### Dealing with temporary positions

E.g. those taken when GPS has powered down so we need to get seasonal/annual displacement from one short survey.

There are likely to be multiple temporary positions in one raw leica file, corresponding to a number of rover locations.

It's probably easiest to do this in a completely separate gps processing folder to avoid filename clashes. Also copy over the relevant base station and orbit files.

First convert the leica file to rinex file(s) using process_rinex. Then window this rinex file into separate rinex files for each rover. Flight timings are very helpful here. You can first check the contents of the rinex file, for example using teqc (only for RINEX2 files):

```bash
teqc +qc rinex_file_name
```

With RINEX3 files, use `gfzrnx` instead:
```bash
gfzrnx -finp <filename> -stk_epo <time_bin_in_seconds>
```
A good time bin is 1800 seconds, i.e. 30 minutes.
  
Check satellite status for each site - if a site is only seeing 4 satellites for a period, or there are significant TRACKing problems, consider windowing out the poor data.  

To do the windowing:

```bash
# Either teqc (RINEX2 only)
teqc -st hhmm00 +dm mm input_rinex_file > output_rinex_file
# Or gfzrnx (RINEX2,3)
gfzrnx -finp <filename> -epo_beg <begin> -d <seconds>
```
  
N.b. in TEQC, option -st by default assumes everything is in seconds, hence you have to specify the number of seconds to force it to understand hours and minutes. So e.g. to extract Lev6 temporary position from 2012 autumn file, begining at 12pm and continuing for 28 minutes:

```bash
teqc -st 120000 +dm 28 levr2430.12o > lev62430.12o
# Equivalent with gfzrnx:
gfzrnx -finp levr2430.12o -epo_beg 2012243_120000 -d 1680 -fout lev62430.12o
```

(Use the same filename format as proper rinex files so that TRACK knows what to look for.)
 
Next, run a PPP analysis on the RINEX file using:

https://webapp.csrs-scrs.nrcan-rncan.gc.ca/geod/tools-outils/ppp.php

You now have two different options:

* Use the PPP analysis to provide APR coordinates to TRACK, then process the RINEX file kinematically with `process_dgps`.
* Use the outputs of the PPP analysis directly, by converting them to a GEOD parquet file of the same format produced by `concatenate_geod`. To do this, adapt the script `ppp_csv_to_parquet.py`.

Copy the GEOD results files into your main GPS processing directory. You'll then be able to concatenate the output GEOD file to the main rover dataset when you run concatenate_geod. 

#### Level-0 data

The RINEX files constitute Level-0 data for archival. Compress these files to CompactRINEX using RNXCRZ, and then gzip them to compress further.


### Concatenate the TRACK files

The GEOD TRACK output files need to be combined together, removing the overlapping hours.

Use `conc_daily_geod` to produce multi-day Parquet files. Each file corresponds to a 'batch' (see the explanation earlier in this readme). These Parquet files are created from data stored in $GNSS_PATH_TRACK_OUT, which is also where the Parquet files get saved to.

In the simple case of collecting data from the field in spring each year and processing using the same base-rover combination, there will be two data batches:

1. rover_base_year1_startDOY_endDOY
2. rover_base_year2_startDOY_endDOY

However, if a different base station needs to be used for certain periods, there will be multiple batches, e.g.

1. rover_bas1_year1_100_200
2. rover_bas2_year1_200_365
3. rover_bas1_year2_1_100


### Export to Level-1 product

We consider the Level-1 data product to be a continuous Parquet file containing all data from all batches since a site was first established. You can generate the Level-1 Parquet file either:

1. Using all the individual data batch Parquet files, or
2. Providing an existing Level-1 Parquet file and only the new data batches.

Use the command line tool `export_level1` in both cases. The Level-1 data are saved in `$GNSS_L1DIR/<site>`.


  
## Detailed Usage: Processing from Level-1 to Level-2: To displacements and velocities

### Installation

As well as installing the package (see very top of this README), you need the following directory structure. 

* `GNSS_L1DIR` : location to store Level-1 data
* `GNSS_L2DIR` : location to store Level-2 data
 	* `{site}`
		* `rotation_{site}.dat` (if already created)
		* `origin_{site}.csv` (if already created)
		* `exclusions_{site}.csv` (if needed)

Set the paths to these folders using environment variables, which are used by the processing scripts to find and export data:

* For Unix environments (including MacOS): see and update the accompanying shell script `gnss.sh`.  You can either transfer these variables into your `.bashrc` (or equivalent), or `source` the shell script.
* For Windows environments: you will need to set equivalent environment variables using the Windows settings. See e.g. https://phoenixnap.com/kb/windows-set-environment-variable for more information.


### Site origin

If older data for this site has already been post-processed then a file, `origin_<site>.csv` will have been generated. Place this file in $GNSS_L2DIR as indicated above.

If this is the first occasion of processing for this site,  run 

```bash
calculate_local_origin <site> <Parquet file>
```

### Correcting pole leans

This can be attempted by fitting functions to remove the lean and leave the residual.
It's best to write your own script to do this on a case by case basis.


### Displacements and velocities of batch(es)

**Level-2 data generation:**  Use `gnss_disp_vel.py`. This can be run on the command line or as a Notebook.

If older data for this site have already been post-processed then a file `rotation_<site>.dat` will exist, defining the coefficients to rotate the coordinates into along/across-TRACK displacement. Place this file in `$GNSS_L2DIR` as indicated above.

Some hints on a workable processing strategy:

* Use the script "iteratively" to identify periods which should be excluded.
* Add exclusion periods to a file named `exclusions_<site>.csv`, located in `$GNSS_L2DIR` (or specify elsewhere with the `-optpath` option of `gnss_disp_vel.py`).
* Re-run the script.
* If corrections due to pole re-drilling are needed then this functionality will first need to be implemented, as it was not (yet) needed for the high-elevation ice velocities campaign!

This script can be used in at least two or three ways:

1. Processing a batch of continuous occupation data. Supply just one parquet file.
2. For a very short data batch - i.e. where a site has been re-occupied for only minutes to hours - use `gnss_disp_vel.py` with the `-stake` option. This disables the smoothing procedures, as they are only applicable to longer time series data. In this case only `xyz` data will be saved to disk.
3. Or let the script process both continuous and daily occupation data automatically. Do not provide `-stake`. Check the messages to make sure that each period has been identified correctly.

Running option 3 within an ipython terminal:

```python
>>> %run gnss_disp_vel.py <site> -f <file1.parquet> <file2.parquet> ...
```

The script applies different filtering and averaging approaches depending on whether a data period has been occupied continuously or only for a short period (e.g. an hour).


### Estimating seasonal and annual displacements

Use `seasonal_annual_disp.py` to produce CSV files per-site of seasonal and annual velocities and displacements.


### Internal archival/retention

Archive the raw receiver files internally. 

Retain the following files generated during processing for internal (re-)use:

- Parquets of GEOD batches, i.e. those output by `conc_daily_geod.py`
- Rotation file: `rotation_<site>.dat`
- Origin file: `origin_<site>.csv`
- Exclusions file: `exclusions_<site>.csv`
- Corrections file: `corrections_<site>.csv` (not currently implemented 2022-07).

### Public data deposition

* Level-0: CompactRINEX files.
* Level-1: Parquet files containing all GEOD data; processing logs.
* Level-2 HDF-5 files, CSV files where relevant; also PNGs.
* Also deposit relevant metadata (installation/maintenance logs, photos).

	
## Other information

### Main differences between Leica and GVT systems

Leica:

- GPS only, therefore IGS FINAL sp3 sufficient
- RINEX2 sufficient, therefore TEQC still works

GVT:

- Multi-constellation, so need multi-c sp3 files
- RINEX3 files, so need to use `gfzrnx`
- Need to tell TRACK which signals to use (with `obs_file`)
- Need to tell TRACK which constellations to use (with `tr_gnss`)


### Telling TRACK which observables to use

If TRACK exits with `**DISASTER** No overlapping data at site` yet the RINEX files are definitely for the same day, then the issue is likely with specifying the observables. See email AT-Mike Floyd, 14.07.2025:

You likely just need to give track a little guidance on what observables to use, since RINEX 3 has so many to choose from.

You will see in your log file that, without any guidance, track has selected the following observables to use:

```
Observables used: 

Site S    Primary      |    Secondary
klsq G L1C L2W C1C C2W |  --- L5Q --- C5Q  
ilhe G L1C --- C1C --- |  --- --- --- ---  
```

If you are only using GPS, [CL]1C and [CL]2W are indeed the most consistently available observables for KLSQ, with [CL]2L being your next most consistently available observables followed by [CL]5Q. But ILHE, using the u-blox receiver, is recording only [CL]2L for its lower frequency, which is the L2C broadcast signal that low-cost receivers are able to record.

So I would simply try adding more information to the "obs_file" option at the top of your command file, e.g. see the track help page (~/gg/help/track.hlp) for details.

Specifically, I would change your current "obs_file" option to something like

```
 obs_file
  <S07> <S07>/<S07>_<day>0.<S09>o F L1C L2W C1C C2W L1C L2L C1C C2L
  <S08> <S08>/<S08>_<day>0.<S09>o K L1C L2L C1C C2L
```

This at least got past the `**DISASTER** No overlapping data at site
ilhe` message and produced a complete result.


### Additional notes about TRACK GEOD files

Column descriptions:

* `# (DD)` becomes `N` when read by the functions in the post-processing module. It specifies the number of double differences per epoch. Per TRACK docs, this should be greater than 0; per Doyle et al. (2014) this could be thresholded at > 4. 
* `NotF` means 'Not Fixed', i.e. the number of ambiguities that have not been fixed to integers in that epoch.
* `#BF` means the number of ambiguities that are needed for the data being analysed.


### Trimble files

Net RS files are in T00 format. You'll need to runpkr00 utility to convert them to dat files first. gnss.RinexConvert.window_overlap can then process the dat file, e.g.:

```python
import gnss
rx = gnss.RinexConvert()
rx.window_overlap("file.dat","22:00:00",28)
```


### Problems with noisy pseudo-range data (lev7)

Email exchange between AJT and MAK, December 2013-January 2014.
"The other thing that can cause this sort of numerical instability is that the site motion is too tightly constrained. It turns out that is what it looks like was resulting in the NaNs. After fiddling for a couple of days it seems likely there are some other issues as well, and I think the receiver is suspect. I've copied below a .cmd file that seemed to produce sensible results. I think the main problem here is that the pseudorange data are very noisy, and as a result too many cycle slips are flagged and incorrect ambiguities are being fixed. To overcome this I changed the noise for the pseudoranges to 10m and modified some thresholds for detecting slips and fixing them to integers. "



## Credits

* Andrew Tedstone 2012-2014 and 2022-25.
* Developed from previous processing routines written by Matt King, Andrew Sole, Ian Bartholomew.

