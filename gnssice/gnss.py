#!/usr/bin/python
"""Kinematic GNSS Processing

Classes:
RinexConvert -- functionality for converting raw leica mdb files to the 
research group rinex file specification. Wrapper to TEQC.
Kinematic -- wrapper to kinematic GPS processing using Track.

Functions:
shellcmd -- used to run subprocesses through shell.
neighborhood -- search behind/current/ahead in list simulataneously (generator).
confirm -- Prompt for yes/no response from user.

Pre-requisites:
As well as all required module imports, TEQC/Gamit/Track must be available
in the environment.

More information:
See the README.

HISTORY
Created on Thu Feb 02 10:50:25 2012 
2022-04: Upgrade to Py3. Some refactoring to Pandas.
2025-07: Major upgrades to support Cryologger GVT/RINEX3.

@author: Andrew Tedstone (andrew.tedstone@unil.ch)

"""
from __future__ import annotations
import matplotlib.pyplot as plt

import numpy as np
import scipy.stats
import subprocess
import datetime
import math
import logging
import os
import pandas as pd
import re
import click
from pathlib import Path
from glob import glob
import shutil

INSTITUTION = 'PLACEHOLDER'
OBSERVER = 'PLACEHOLDER'
GVT_FW = 'v1.32' # GVT firmware version, to put in RINEX headers

# Look-up table taken from sh_get_orbits
ORBIT_PRODUCTS_TO_SP3_NAME = {
   'igsf':'igs',
   'igsr':'igr',
   'igsu':'igu',
   'codf':'cof',
   'code':'cod',
   'codm':'com',
   'codr':'cor',
   'emrf':'emr',
   'esaf':'esa',
   'gfzf':'gfz',
   'gfzm':'gbm',
   'grgm':'grm',
   'jaxm':'jam',
   'jplf':'jpl',
   'mitf':'mit',
   'mitm':'mim',
   'ngsf':'ngs',
   'siof':'sio',
   'sior':'sir',
   'siou':'siu',
   'tumm':'tum',
   'wuhm':'wum'
}

# Look-up table of receivers to TRACK command filename labels
# Key: the name of the receiver encoded in the RINEX header block
# Value: the rcvr ID to substitute into track_<base>_<rcvr>.cmd.
TRACK_CMD_RCVR = {
    'LEICA GX1220+':'L1200',
    'CryoLogger GVT':'gvt'
}

def shellcmd(cmd, timeout_seconds=False, retry_n=2):
    """A general wrapper to the Popen command, running through the shell.
    
    Parameters
    ----------
    cmd : the shell command to run.
    timeout_seconds : Make the process timeout after n seconds. Optional. 
    retry_n : number of times to retry the process if it times out. Only valid 
    if timeout_seconds is set. Optional.
    
    Returns
    -------
    out : dictionary with keys 'stdout' and 'stderr'. 
    
    """
    if timeout_seconds != False:
        retry_count = 0
        retry = True
        while retry == True and retry_count < retry_n:
            retry_count += 1
            try:
                pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
                stdout,stderr = pid.communicate(timeout=timeout_seconds)                    
            except subprocess.TimeoutExpired:
                if retry_count < retry_n:
                    retry = True
                    print('Process timed out. Retrying...')
                else:
                    stdout = "Process timed out."
                    stderr = "Process timed out."
                    retry = False
            else:
                retry = False
    else:
        pid = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
        stdout,stderr = pid.communicate()
    
    if type(stdout) is bytes:
        try:
            stdout = stdout.decode('ascii')
        except UnicodeDecodeError:
            print('Warning: could not ASCII-decode stdout')
    if type(stderr) is bytes:
        try:
            stderr = stderr.decode('ascii')
        except UnicodeDecodeError:
            print('Warning: could not ASCII-decode stdout')

    toret = {}
    toret['stdout'] = stdout
    toret['stderr'] = stderr
    return toret

        
def neighborhood(iterable):
    """Enables return of previous and next items in list. Used by read_orbits.
    
    From http://stackoverflow.com/questions/323750/how-to-access-previous-next-element-while-for-looping
    """
    iterator = iter(iterable)
    prev = None
    item = iterator.__next__()  # throws StopIteration if empty.
    for next in iterator:
        yield (prev,item,next)
        prev = item
        item = next
    yield (prev,item,None)


def confirm(prompt=None, resp=False):
    """prompts for yes or no response from the user. Returns True for yes and
    False for no.

    'resp' should be set to the default value assumed by the caller when
    user simply types ENTER.

    >>> confirm(prompt='Create Directory?', resp=True)
    Create Directory? [y]|n: 
    True
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y: 
    False
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y: y
    True
    
    From http://code.activestate.com/recipes/541096-prompt-the-user-for-confirmation/

    """
    
    if prompt is None:
        prompt = 'Confirm'

    if resp:
        prompt = '%s [%s]|%s: ' % (prompt, 'y', 'n')
    else:
        prompt = '%s [%s]|%s: ' % (prompt, 'n', 'y')
        
    while True:
        ans = input(prompt)
        if not ans:
            return resp
        if ans not in ['y', 'Y', 'n', 'N']:
            print('please enter y or n.')
            continue
        if ans.lower() == 'y':
            return True
        if ans.lower() == 'n':
            return False    

def gpsweekD(yr, doy, wkday_suff=False):
    """
    Convert year, day-of-year to GPS week format: WWWWD or WWWW

    Verbatim from https://github.com/GeoscienceAustralia/gnssanalysis/gnssanalysis/gn_datetime.py
    Which is based on code from Kristine Larson's gps.py
    https://github.com/kristinemlarson/gnssIR_python/gps.py

    Input:
    yr - year (int)
    doy - day-of-year (int)
    wkday_suff : If True, return WWWWD, otherwise WWWW

    Output:
    GPS Week in WWWWD format - weeks since 7 Jan 1980 + day of week number (str)
    """

    # Set up the date and time variables
    yr = int(yr)
    doy = int(doy)
    d = datetime.strptime(f"{yr}-{doy:03d} 01", "%Y-%j %H")

    wkday = d.weekday() + 1

    if wkday == 7:
        wkday = 0

    mn, dy, hr = d.month, d.day, d.hour

    if mn <= 2:
        yr = yr - 1
        mn = mn + 12

    JD = np.floor(365.25 * yr) + np.floor(30.6001 * (mn + 1)) + dy + hr / 24.0 + 1720981.5
    GPS_wk = int(np.floor((JD - 2444244.5) / 7.0))

    if wkday_suff:
        return str(GPS_wk) + str(wkday)
    else:
        return str(GPS_wk)

@click.command()
@click.argument('year')
@click.argument('start_doy')
@click.argument('end_doy')
@click.option('--orbit', default='codm')
@click.option('--overlap', default=False)
@click.option('--clearup', default=True)
@click.option('--rename', default=False, help='if True, then year signifier will be removed from filename.')
def get_orbits(
    year: int, 
    start_doy: int, 
    end_doy: int, 
    orbit: str='codm',
    download=True,
    overlap: bool=False,
    clearup: bool=True
    ) -> None:
    """Download daily sp3 orbit files and optionally overlap them.
    
    Each resulting overlapped file contains the previous day, the current
    day, and the next day. (Except for 1 Jan or 31 Dec).
    
    Do not try to download orbits spanning two years - only run on a 
    per-yearly basis.
    
    If clearup=True, the un-overlapped files will be deleted afterwards.

    Files get renamed to OOOYYDDD0.sp3, where OOO=orbit SP3 name, YY=year, DDD=julian day of year.

    """
    start_doy = int(start_doy)
    end_doy = int(end_doy)
    n_days = end_doy - start_doy + 1 
    if n_days < 0:
        print("It looks like you've entered the number of days to download, not the end day. You need to specify the end day. Exiting...")
        return
    Path(os.environ['GNSS_PATH_SP3_DAILY']).mkdir(exist_ok=True)

    def make_target_fn(orbit, year, doy):
        """ Generate the filename format needed for our work flow. """
        return '{o}{y}{d}0.sp3'.format(
            o=ORBIT_PRODUCTS_TO_SP3_NAME[orbit],
            y=str(year)[:-2],
            d=str(doy).zfill(3)
        )

    if download:
        print('Downloading...')
        shellcmd("cd " + os.environ['GNSS_PATH_SP3_DAILY'] )
        cmd = f"sh_get_orbits -orbit {orbit} -yr " + str(year) + " -doy " + \
            str(start_doy) + " -ndays " + str(n_days) + " -nofit"
        print(cmd)
        status = shellcmd(cmd)
        print(status['stdout'])
        print(status['stderr'])
        shellcmd('cd ' + os.environ['GNSS_WORK'])

        print('Renaming...')
        for doy in range(start_doy, end_doy+1):
            old_fn = '{o}{wd}.sp3'.format(
                o=ORBIT_PRODUCTS_TO_SP3_NAME[orbit],
                wd=gpsweekD(year, doy, wkday_suff=True)
            )
            target_fn = make_target_fn(orbit, year, doy)
            os.rename(
                os.path.join(os.environ['GNSS_PATH_SP3_DAILY'], old_fn),
                os.path.join(os.environ['GNSS_PATH_SP3_DAILY'], target_fn) 
            )
    
    if overlap:
        print('Overlapping...')
        Path(os.environ['GNSS_PATH_SP3_OVERLAP']).mkdir(exist_ok=True)
    
        for doy in range(start_doy, end_doy+1):
            tday_fn = os.path.join(os.environ['GNSS_PATH_SP3_DAILY'], make_target_fn(orbit, year, doy))
            yday_fn = os.path.join(os.environ['GNSS_PATH_SP3_DAILY'], make_target_fn(orbit, year, doy-1))
            tomo_fn = os.path.join(os.environ['GNSS_PATH_SP3_DAILY'], make_target_fn(orbit, year, doy+1))

            out_fn = os.path.join(os.environ['GNSS_PATH_SP3_OVERLAP'], make_target_fn(orbit, year, doy))

            cmd = 'cat {prev} {today} {tomo} > {fn}'.format(
                prev=yday_fn,
                today=tday_fn,
                tomo=tomo_fn,
                fn=out_fn
            )
            sout, serr = shellcmd(cmd)
            if serr is not None:
                raise IOError(serr)
        
    if clearup:
        print("Clearing up...")
        shellcmd("rm {1}/*.sp3".format(os.environ['GNSS_PATH_SP3_DAILY']))
    
    print("Done.")


def read_track_geod_file(
    fname: str
    ) -> pd.DataFrame:
    """ 
    Read a GEOD file created by TRACK.

    Note that the last column, 'Flag', is dropped.

    """
    # Specify the header manually - the file contains column names with 
    # white space which Pandas interprets as signifying a new column.
    # Note that some columns are renamed slightly to remove non-alphabet characters.
    names = ['YY', 'DOY', 'Seconds', 'Latitude', 'Longitude', 'Height', 
        'SigN', 'SigE', 'SigH', 'RMS', 'N', 'Atm', 'plus_minus', 'Fract_DOY', 
        'Epoch', 'NBF', 'NotF', 'Flag']

    bounds = [
        (0, 5),
        (7, 10),
        (13, 25),
        (28, 40),
        (42, 55),
        (58, 67),
        (68, 73),
        (74, 79),
        (80, 85),
        (86, 91),
        (92, 94),
        (95, 103),
        (104, 112),
        (115, 130),
        (131, 137),
        (138, 141),
        (142, 145),
        (146, 148)
    ]

    data = pd.read_fwf(
        fname, 
        skiprows=[0, 1],
        na_values=['*'],
        names=names,
        colspecs=bounds,
    )

    data = data.drop(labels=['Flag'], axis=1)
    data = data.apply(pd.to_numeric)
    return data


def read_track_neu_file(
    fname: str
    ) -> pd.DataFrame:
    """ 
    Read a NEU file created by TRACK.

    """
    # Specify the header manually - the file contains column names with 
    # white space which Pandas interprets as signifying a new column.
    # Note that some columns are renamed slightly to remove non-alphabet characters.
    names = ['YY', 'MM', 'DD', 'HR', 'MIN', 'Sec', 'dNorth', 'dNorth_plus_minus',
            'dEast', 'dEast_plus_minus', 'dHeight', 'dHeight_plus_minus', 'RMS',
            'N', 'Atm', 'Atm_plus_minus', 'Fract_DOY', 'Epoch', 'BF', 'NotF', 
            'Flag', 'Rho_UA']

    bounds = [
        (0, 5),
        (6, 8),
        (9, 11),
        (12, 14),
        (15, 17),
        (19, 28),
        (29, 43),
        (44, 53),
        (54, 69),
        (70, 79),
        (80, 95),
        (96, 105),
        (106, 114),
        (115, 118),
        (119, 128),
        (129, 137),
        (140, 155),
        (156, 162),
        (163, 166),
        (167, 170),
        (171, 172),
        (173, 189)
    ]

    data = pd.read_fwf(
        fname, 
        skiprows=[0, 1],
        na_values=['*'],
        names=names,
        colspecs=bounds,
    )

    data = data.drop(labels=['Flag'], axis=1)
    data = data.apply(pd.to_numeric)
    return data
    

class RinexConvert:
    """Accomplish Rinex file processing from binary receiver files."""
    
    def __init__(self):
        """Initialise class variables."""
        self.institution = INSTITUTION
        self.observer = OBSERVER

    def parse_rinex2_filename(fn, to_numeric=False):
        """ Extract and return elements of a RINEXv2 file name 

        to_numeric : bool. If True, parse DOY as integer.
        returns : dict {site, DOY, yr (2-charac), year (YYYY)}
        """
        fn = os.path.split('/')[-1]
        site = fn[0:4]
        doy = fn[4:7]
        yr = fn[-4:-2]
        if to_numeric:
            doy = int(doy)
        year = int(yr) + 2000
        return {
            'site':site,
            'DOY':doy,
            'yr': yr,
            'year':year
        }
        
    def gvt_to_rinex(
        self, 
        input_file : str, 
        site : str,
        rcvr : str=f'NA / CryoLogger GVT / ZED-F9P FW{GVT_FW}'
        antenna : str='NA / SFETOP106',
        site_type : str='Glacier',
        observer : str=None,
        ) -> None:
        """ Convert a daily GVT u-blox file to a daily RINEX3 file. Wraps convbin.

        input_file : the path to the .ubx file.
        site : the 4-character site identifier.
        rcvr : sets the receiver in the RINEX header.
        antenna : sets the antenna in the RINEX header.
        site_type: sets the marker type in the RINEX header.
        observer : sets the observer in the RINEX header.

        Behaviour:
        - We produce RINEX files with 10-second sampling.
        - We allow a tolerance of 0.01 seconds.
        - We use the receiver-specific option TADJ to force floating-point seconds to integer values.
        - We output filenames in the RINEX2 format, for backward compatibility with existing workflow.

        """
        if observer is None:
            observer = '{o} / {i}'.format(i=self.institution, o=self.observer)

        output_dir = os.path.join(os.environ['GNSS_PATH_RINEX_DAILY'], site)
        Path(output_dir).mkdir(exist_ok=True)

        cmd = ( f'convbin -hm {site} -c {site} -ho {observer} -hr {rcvr} '
                f'-ha {antenna} -ht {site_type} -ti 10 -tt 0.01 -ro "TADJ=1.0 '
                f'-o \%r\%n0.\%yo -d {output_dir} {input_file}'
            )
        stdout, stderr = shellcmd(cmd)

        return (stdout, stderr)


    def leica2rinex(
        self, 
        input_file : str, 
        site : str,
        output_file : str
        ) -> None:
        """Convert a leica mdb file to a rinex file. A simple wrapper to teqc.
        
        Whilst this will convert very large single mdb files, be aware
        that the resulting rinex files may not be usable! There appears to be
        a size limit, maybe around 2.5Gb. You can use window_overlap to produce
        daily rinex files from the single mdb file instead.
        """
        cmd = "teqc -leica mdb -O.o '" + self.observer + "' -O.ag '" + \
        self.institution + "' " + input_file + " > " + output_file
        status = shellcmd(cmd)
        print(status)
                
    
    def leica_joindaily(self, input_file_prefix, site,moving=1):
        """Convert daily leica files to rinex and join together into one big file.
        
        If the series of files contain data from more than one year, one big
        file will be generated for each file, according to the rinex filenaming
        convention.
        
        The series of files (i.e. all those with the same file_prefix) must all
        have been recorded at the same site (monument). The site name will
        be set as the monument name in the resulting rinex files.
            
        Inputs:
            input_file_prefix: the filename common to all files in the series - it
            comes before the file ending, i.e. <file_prefix>.YYo. Required.
            site: 4-letter site identifier. Required. 
            moving: 1 if site moving, 0 if static. Default is 1.
        Outputs:
            joined_files: list of filenames of created joined files.
        """
        # Record the years being examined in this file batch        
        years = []
        # Record any files which teqc has problems with
        problems = []
        # Get list of daily files to convert
        filestring = shellcmd("ls " + input_file_prefix + ".*")
        filestring = filestring['stdout']
        files = st.split(filestring, sep='\n')
        for fn in files:
            if fn == '':
                continue
            # Extract file date and convert to doy
            # Rinex returns file start datein format: 
                # '-O.st[art] yyyy mm dd hh mm s.sss...'
            rinex_st = shellcmd("teqc ++config " + fn + r" | grep '\-O\.st\[art\]' -a")
            if rinex_st['stderr'] != '':
                if rinex_st['stderr'].find("Notice") < 0:
                    print("Couldn't grep start from ++config: " + rinex_st['stderr'] + ', continuing...')
                    problems.append(fn)
                    continue
            rinex_st = rinex_st['stdout']
            doy = datetime.datetime(int(rinex_st[11:15]),int(rinex_st[16:18]),
                                    int(rinex_st[19:21])).timetuple().tm_yday
            print(fn + " (day " + str(doy) + ")")
            # Update list of years to make large files from
            if int(rinex_st[13:15]) not in years:
                years.append(int(rinex_st[13:15]))
            # Run teqc            
            output_fn = site + "_" + str(doy) + "0_d." + rinex_st[13:15] + "o"        
            cmd = "teqc -leica mdb -O.mo " + site + " -O.o '" + self.observer + \
            "' -O.ag '" + self.institution + "' -O.mov " + str(moving) + " " + fn + " > " + output_fn
            status = shellcmd(cmd)
            if status['stderr'] != '' and status['stderr'].find("Notice") == -1:
                print("Teqc conversion failed: " + status['stderr'] + ', continuing...')
                problems.append(fn)
            else:
                print ("   ...processed")            
            
        # Combine all files into one for each year
        print("Starting to join rinex files together. This might take a while.")
        joined_files = []
        for y in years:
            fn = "all_" + site + "." + str(y) + "o"
            cmd = "teqc " + site + "_*." + str(y) + "o > " + fn
            status = shellcmd(cmd)
            print(status['stdout'])
            print(status['stderr'])
            joined_files.append(fn)
            
        if len(problems) > 0:
            print("There were problems processing some files, \
            check list of returned filenames.")
       
        toret = {}
        toret['problems'] = problems
        toret['joined_rinex'] = joined_files
        return toret
        

    def leica_to_daily_rinex(
        self, 
        input_file : str,
        site: str=None, 
        moving: int=1,
        start_date: datetime.datetime=None,
        end_date: datetime.dateime=None
        ) -> None:
        """ 
        Generate daily RINEX files from one large Leica MDB file.
        
        Inputs:
            input_file: filename of the rinex file to split into windowed files.
            site: recording site 4-character identifier.  
            moving: Is the site moving or static? 1 if moving, 0 if static. 
                    Defaults to 1 (moving)
            start_date: date to start processing from. 
            end_date: date to end processing at. 
            Supply dates as datetimes.
            
        E.g. Process from a raw leica file:
            leica_to_daily_rinex("leica-file.m00","0",24,
                                 site="levb", 
                                 start_date=dt.datetime(...),end_date=dt.datetime(...))

        Outputs:
            None.
            
        """

        add_leica = "-leica mdb -O.o '" + self.observer + "' -O.ag '" + \
        self.institution + "' -O.mo '" + site + "' -O.mov " + str(moving) + " "

        # This returns today's Julian day of Year ('yday' == 'year day', not yesterday!)
        start_doy = start_date.timetuple().tm_yday
        end_doy = end_date.timetuple().tm_yday
        
        print("Commencing windowing on " + input_file)
        print("This file begins on " + start_date.strftime("%Y.%m.%d (DOY %j)"))
        print("and ends on " + end_date.strftime("%Y.%m.%d (DOY %j)"))

        if st_offset != '00:00:00':
            # The window commences the day before the main day of interest
            cal_date = start_date - datetime.timedelta(days=1)
        else:
            # Window will start on day of interest
            cal_date = start_date
        
        out_root = os.path.join(os.environ['GNSS_PATH_RINEX_DAILY'], site)
        Path(out_root).mkdir(exist_ok=True)

        for doy in range(start_doy,end_doy + 1):
            print("Processing day " + str(doy))
            output_file = os.path.join(
                out_root,
                '{site}{doy}0.{yr}o'.format(site=site, doy=str(doy).zfill(3), yr=cal_date.strftime("%y"))
                )
            cmd = "teqc " + add_leica + "-st " + cal_date.strftime("%Y%m%d") + " +dh 24" + \
            " " + input_file + " > " + output_file
            print("    " + cmd)
            status = shellcmd(cmd)
            if status['stderr'] != '': #and status['stderr'].find("Notice") == False:
                print("         " + status['stderr'])
            # Increment date.
            cal_date = cal_date + datetime.timedelta(days=1)

        print("Done.")


    def window_overlap(
        self, 
        input_file : str,
        st_timestart : str='000000', 
        dh : int=24
        ):
        """ Make daily RINEX file with overlapping window into the preceding and subsequent days.

        input_file : str with RINEX2 filename. This function assumes that the DOY-1 and DOY+1 files are located in the same place.
        st_timestart : time at which to begin the overlap, following syntax allowed by gfzrnx.
        dh : window length in hours

        """

        # Get info about this file / date
        finfo = self.parse_rinex2_filename(input_file, to_numeric=True)
        doy_prev = str(finfo['DOY'] - 1).zfill(3)
        doy_next = str(finfo['DOY'] + 1).zfill(3)
        doy = str(finfo['DOY']).zfill(3)
        yr = finfo['yr']

        # Construct the command line parameters
        f_prev = os.path.join(, args.site, f'{site}{doy_prev}0.{yr}o')
        f_next = os.path.join(os.environ['GNSS_PATH_RINEX_DAILY'], args.site, f'{site}{doy_next}0.{yr}o')
        f_out = os.path.join(os.environ['DATA_PATH_RINEX_OVERLAP'], args.site, f'{site}{doy}0.{yr}o')
        epo_beg = str(finfo['year']) + doy + '_' + st_timestart
        seconds = str(60 * 60 * dh)

        # Construct the command
        cmd = ( 
            f'gfzrnx -finp  {f_prev} {f} {f_next} '
            f'-fout {f_out} '
            f'-epo_beg {epo_beg} -d {seconds}'
        )
        sout, serr = shellcmd(cmd)     
        return(sout, serr)
        
        
       
class Kinematic:
    """Kinematic GNSS processing functionality, utilising Track and Gamit."""
    
    def __init__(self):
        self.ion_stats = None
        self.MW_WL = None
        self.LG = None
        self.apriori = None
        self.config_subfolder = ""
         

    def track(self, base, rover, 
        year, doy_start, doy_end, 
        show_plot=True, 
        use_auto_qa=True, spearman_threshold=None, rms_threshold=22.5):
        """Wrapper to track kinematic processing.
        
        Processing takes one of two slightly different approaches. With 
        use_auto_qa=True, each day of data will be checked according to one of two
        methods:
        (1) linearity in East versus North using Spearman correlation. If the result of the test
        is > spearman_threshold, the day is accepted automatically. 
        (2) Median RMS from TRACK below rms_threshold.
         
        If the test passes then processing moves on to the next day with no interaction necessary. 
        If the day fails the test, the user must choose the next step to take.
        
        If use_auto_qa=False, the user must approve every day processed manually.
        
        The function must be used interactively with either approach - it does
        not work in a batch environment.
        
        A log file for each rover is also created which contains a summary of the
        processing parameters used, together with user assessment of quality
        and any further comments. 'A' and the Spearman test coefficient are 
        output if use_auto_qa=True and the day passes the test.
        
        Note that any -ve signs on APR coordinates extracted from the previous
        day's files will be removed as track cannot cope with these being
        provided in the command line arguments. The negative sign must
        be put in the site-specific track cmd file.
        
        Cmd files are named as track_<base>.cmd
        
        Inputs:
            base: 4-character identifier of static base station.
            rover: 4-character identifier of moving (on-ice?) GPS station.
            year: year of data to process.
            doy_start: day of year to start on.
            doy_end: day of year to end on.
            show_plot: if True, plot of data will be popped up automatically.
            use_auto_qa: if True, high quality days will automatically be 
              accepted. When a day is automatically accepted, a plot of the 
              data will be saved but will not be displayed on screen, even 
              if show_plot=True.
           spearman_threshold: the value that the spearman coefficient must 
              exceed for the day to be approved automatically. Pass None to disable.
            rms_threshold: value in mm that median RMS must not exceed for day to
            be approved automatically. Pass None to disable.
        Outputs:
            None.
            
        """

        # Check we are in GNSS working directory
        if os.getcwd() != os.environ['GNSS_WORK']:
            raise RuntimeError('Your current working directory does not match GNSS_WORK.')
                
        # Set up logging. 
        lfn = "gps.track." + rover + ".log"
        logging.basicConfig(filename=lfn,level=logging.INFO,format='%(message)s')
        logging.info("\n\n" + datetime.datetime.now().strftime("%Y-%m-%d %I:%M:%S %p") +
        ": BEGINNING NEW PROCESSING BATCH")
        logging.info("Year: " + str(year) + ". Start: day " + str(doy_start) + ", End: day " + \
        str(doy_end) + ", Base/Static: " + base)

        yr_short = year % 1000

        # Identify rover receiver type from RINEX
        fpth = os.path.join(os.environ['GNSS_WORK'], 'rinex', rover, f'{rover}{doy_start}0.{yr_short}o')
        rcv_type, serr = shellcmd(f"grep -E 'REC # / TYPE / VERS' {fpth} | awk '{print $2, $3}'")
        rcvr = TRACK_CMD_RCVR[rcv_type]

        # Now build TRACK cmd filename
        cmdfile = os.path.join(os.environ['GNSS_WORK'], self.config_subfolder, f'track_{base}_{rover}_{rcvr}.cmd')
        # Check it is exists before we go any further
        if not os.path.exists(cmdfile):
            raise IOError(f'Track CMD file {cmdf} not found.')

        # Set up output directory
        output_dir = os.path.join(os.environ['GNSS_PATH_TRACK_OUT'], rover)
        Path(output_dir).mkdir(exist_ok=True)

        # Enter main processing loop, works on a per-day basis
        for doy in range(doy_start,doy_end):

            save_opts = dict(r=rover, b=base, y=year, d=str(doy).zfill(3))

            print("Processing day " + str(doy).zfill(3) + "...")
            # Deal with APR coordinates            
            if doy == doy_start and self.apriori != None:
                print("Using specified APR coordinates.")
                apriori = str(self.apriori[0]) + " " + str(self.apriori[1]) + \
                " " + str(self.apriori[2])
                logging.info("Using user-input APR coordinates: " + apriori)
            else:
                print("Extracting APR coordinates from yesterday's data.")
                try:                
                    lc = open("track.GEOD." + rover + ".LC","r")
                except IOError:
                    print("!!APR .LC file does not exist! Terminating processing.")
                    logging.info("Day " + str(doy).zfill(3) + ": APR .LC file does not exist, terminating")
                    return False                    
                lines = lc.readlines()
                lc.close()
                last_line = lines[-1]
                vals = last_line.split()
                # Convert geodetic to cartesian
                # File columns are in order 3=lat, 4=lon, 5=height
                # Note that convertc takes order lon lat height.
                status = shellcmd("convertc " + str(vals[4]) + " " + 
                str(vals[3]) + " " + str(vals[5]) + " XYZ")                
                xyz = status['stdout'].split()
                apriori = str(xyz[0]).strip('-') + " " + str(xyz[1]).strip('-') + " " + str(xyz[2]).strip('-')
            
            # This loop enables reprocessing with changed parameters
            retry = False
            while True:
                ion_stats = self.ion_stats    
                MW_WL = self.MW_WL
                LG = self.LG
                exclude_svs = ''                
                if retry == True:
                    print("Reprocessing day...enter new values or press Return to use Default.")
                    ion_stats = input("    @ion_stats <jump>: ")
                    if len(ion_stats) == 0: ion_stats = self.ion_stats
                    MW_WL = input("    @float_type <WL_Fact> (MW_WL weighting): ")
                    if len(MW_WL) == 0: MW_WL = self.MW_WL
                    LG = input("    @float_type <Ion_fact> (LG Combination Weighting): ")
                    if len(LG) == 0: LG = self.LG
                    print("Processing with new parameter values...")
                else:
                    print("Processing with defaults...")
                        
                # Construct track argument.
                # Parameters (-s) are in exact order expected by the cmd file, 
                # do not change!!
                outf = rover + "_" + base + "_" + str(year) + "_" + str(doy).zfill(3) + ".out"
                
                cmd = 'track -f {cf} -d {d} -s {ap} {mw} {lg} {istat} {base} {rover} {yr} > {out}'.format(
                    cf=cmdfile,
                    d=save_opts['d'],
                    ap=apriori,
                    mw=MW_WL,
                    lg=LG,
                    istat=ion_stats,
                    base=base,
                    rover=rover,
                    yr=yr_short, 
                    out=outf
                    )
                print(cmd)

                # Send to track. 
                status = shellcmd(cmd,timeout_seconds=1200,retry_n=1)
                
                # Check track status, this catches non-IOSTAT errors (e.g. SP3 Interpolation errors)
                if status['stderr'] != '':
                    plt.title('ERROR - track terminated. Close this figure window and respond to command prompt.')
                    plt.show()
                    print("ERROR: Track terminated with the following error message:")
                    print(status['stderr'])
                    track_error = True
                    while True:
                        action = input("[T]ry again, [S]kip day, [H]alt processing session?: ").upper()
                        if action in ['T','S','H']:
                            break
                        else:
                            print('Not a valid option.')
                    if action == 'T': # Try again
                        logging.info("!Day " + str(doy).zfill(3) + ": Track failed to process. Retrying day. Error: " + str(status['stderr']))
                        retry = True
                        # Go through another iteration of the reprocessing loop
                        continue 
                    elif action == 'S': # Skip day
                        print('Skipping day...')
                        logging.info("!Day " + str(doy).zfill(3) + ": Track failed to process. Skipping day. Error: " + str(status['stderr']))
                        retry = False
                        # Break out of the reprocessing loop
                        break
                    elif action == 'H': # Halt processing session
                        print('Processing halted.')
                        exit()
                    
                                                            
                # Check for IO errors, otherwise show RMS values.
                fid = open(outf)
                track_error = False
                store_rms = []
                for line in fid.readlines():
                    if "IOSTAT error" in line:
                        print("Track IOSTAT error: " + line)
                        print("...skipping this day.")
                        logging.info("!Day " + str(doy).zfill(3) + ": Track failed to process.")
                        track_error = True                        
                        # Break out of this loop                        
                        break
                    if "Average RMS" in line:
                        print(line.strip("\n"))
                        store_rms.append(float(re.search(r'[0-9]+\.[0-9]+', line).group()))
                fid.close()
                if track_error == True:
                    break # break out of while:true loop
                
                # Do automated quality check, if requested.
                if use_auto_qa == True: 
                    if spearman_threshold != None:
                        data = read_track_neu_file("track.NEU." + rover + ".LC")                      
                        spearman = scipy.stats.spearmanr(data['dEast'], data['dNorth'])           
                        print('Spearman value: ' + str(spearman[0]))
                        if spearman[0] < 0:
                            spearman_v = spearman[0] * -1
                        else:
                            spearman_v = spearman[0]
                        if spearman_v > spearman_threshold:
                            keep = True
                            show_plot = False
                        else:
                            print('Day rejected by Spearman test.')
                            keep = False
                            show_plot = True
                        comment = 'Spearman: %s' %spearman[0]
                    elif rms_threshold != None:
                        rms = np.median(np.array(store_rms))
                        print('Median RMS: %s' %rms)
                        if rms < rms_threshold:
                            keep = True
                            show_plot = False
                        else:
                            print('Rejected by RMS.')
                            keep = False
                            show_plot = True
                        comment = 'Med.RMS: %s' %rms
                
                # If day passed automatically, prevent the plot from popping up
                if use_auto_qa == True and keep == True:
                    __show_plot = False
                # Otherwise continue to use user-specified option.
                else:
                    __show_plot = show_plot

                # Save a scatter plot, also display subject to above.
                plot_fname = "track.NEU." + rover + ".LC"
                save_plot_to = '{r}_{b}_{y}_{d}_{ftype}.png'.format(ftype='NEU', **save_opts)
                view_track_output(base, rover, doy, 
                                fname=plot_fname,
                                display=__show_plot,
                                save_to=save_plot_to)
                
                # Do manual quality check if automatic not on or if automatic test failed.
                if use_auto_qa == False or (use_auto_qa == True and keep == False):
                    keep = confirm("Keep these results? (press Enter to accept, n to reject ",resp=True)
                    if keep == False:
                        print("Day rejected.")
                    comment = input("Comment for logging (optional): ")
                    if keep == True:
                        quality = input("Quality indication (good=G,ok=O,bad=B): ").upper()
                    else: 
                        quality = "REJECTED"
                
                # Set parameters for logging if automated outputs a keeper
                elif keep == True:
                    quality = 'A'
                
                # Sanity check
                else:
                    print('Why are we here?')
                
                # Prepare and save log entry
                log_str = "Day " + str(doy).zfill(3) + "     ion_stats=" + str(ion_stats) + \
                " MW_WL=" + str(MW_WL) + " LG=" + str(LG) + " | Q:" + \
                str(quality) + " | " + comment
                logging.info(log_str)
                
                if keep == True:
                    # Break out of the while True loop
                    break;
                else:
                    retry = True
            
            if track_error == False:
                print("Saving these results.")
                # NEU
                os.rename(
                    'track.NEU.{1}.LC'.format(rover), 
                    os.path.join(output_dir, '{r}_{b}_{y}_{d}_NEU.dat'.format(**save_opts))
                    )
                # GEOD (we copy rather than move this as the .LC file is still needed for extracting APR coordinates)
                shutil.copy(
                    'track.GEOD.{1}.LC'.format(rover), 
                    os.path.join(output_dir, '{r}_{b}_{y}_{d}_GEOD.dat'.format(**save_opts))
                    )
                # Plot of results
                os.rename(
                    save_plot_to, 
                    os.path.join(output_dir, save_plot_to)
                    )
                # Track OUT file
                os.rename(
                    '{r}_{b}_{y}_{d}.out'.format(**save_opts), 
                    os.path.join(output_dir, '{r}_{b}_{y}_{d}.out'.format(**save_opts))
                    )
            
        print("Batch finished.")
    
    
# @click.command()
# @click.argument('base')
# @click.argument('rover')
# @click.argument('doy')
# @click.option('--gtype', default='NEU', help='currently only NEU supported')
def view_track_output(base, rover, doy, gtype='NEU', fname=None, 
                      display=True, save_to=None):
    """Display a scatter plot of reconciled daily track data.
    
    Inputs:
        base: 4-character identifier of base (static) receiver
        rover: 4-character identifier of moving receiver
        doy: day of year to examine
        gtype: georeferencing type. Optional. Currently only NEU supported
        fname: filename of file. Optional. If not used, the filename will
            be deduced from the other input parameters.
        save_to: None (don't save), or path/filename to save plot to, including extension.
    
    Outputs:
        None.
    
    """
    if fname == None:        
        fname = rover + '_' + base + '_' + str(doy) + gtype + '.dat'
    
    if gtype=='NEU':    
        data = read_track_neu_file(fname)
    else:
        print('Files georeferenced in a format other than NEU are currently unsupported.')
        return
    
    plt.figure()
    
    plt.plot(data['dEast'], data['dNorth'],
                c='tab:blue', marker='.', linestyle='none',
                label='Positions', alpha=0.2)
    only_fixed = data[data['NotF'] == 0]
    plt.plot(only_fixed['dEast'], only_fixed['dNorth'],
                c='tab:orange', marker='.', linestyle='none',
                label='Position (unfixed epoches removed)', alpha=0.2)
    plt.plot(data['dEast'].iloc[0],data['dNorth'].iloc[0],
             c='tab:red', marker='+', markersize=10,
             label='Start')
    plt.plot(data['dEast'].iloc[-1],data['dNorth'].iloc[-1],
             c='tab:red', marker='+', markersize=10,
             label='End')
    
    plt.ylabel('North (m)')
    plt.xlabel('East (m)')
    plt.grid()
    plt.title('%s %s %s' %(base, rover, doy))
    #plt.axes().set_aspect('equal','datalim')    
    plt.legend(numpoints=1, scatterpoints=1, markerscale=0.6, loc='best')    
    if save_to is not None:
        plt.savefig(save_to, orientation='landscape', dpi=200)  
    if display == True:
       plt.show()
    else:
        plt.close()
    return             