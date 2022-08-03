#!/usr/bin/python
"""Kinematic GPS Processing

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

@author: Andrew Tedstone (andrew.tedstone@unifr.ch)

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

def shellcmd(cmd,timeout_seconds=False,retry_n=2):
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
    """Accomplish Rinex file processing from Leica MDB files."""
    
    def __init__(self):
        """Initialise class variables."""
        self.institution = "University of Fribourg"
        self.observer = "Andrew Tedstone"
        

    def leica2rinex(self, input_file, site, output_file):
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
        
                 
    def window_overlap(self, input_file, st_offset, dh,
                       leica=False,
                       site=False, moving=1,
                       start_date=False, end_date=False):
        """ Split one large rinex file into windowed files.
            
        Inputs:
            input_file: filename of the rinex file to split into windowed files.
            st_offset: The time by which to offset from midnight. Include minus
                    sign if necessary. Eg. -22:00:00 starts processing at 10pm 
                    on day before. Provide as a string, i.e. "-22:00:00"
            dh: duration of window in hours. E.g. 28. Integer
            leica: set to True if input file is in raw leica format. If True,
                also provide...
                    site: recording site 4-character identifier.  
                    moving: Is the site moving or static? 1 if moving, 0 if static. 
                    Defaults to 1 (moving)
                    start_date: date to start processing from. 
                    end_date: date to end processing at. 
                    Supply dates in format [yyyy,m,d].
            
        E.g. Process from a raw leica file:
            window_overlap("leica-file.m00","-22:00:00",28,
                                 leica=True,
                                 site="levb", 
                                 start_date=dt.datetime(...),end_date=dt.datetime(...))
        E.g. Process from a rinex file:
            window_overlap("rinexfile.11o","-22:00:00",28)
        
        
        Outputs:
            None.
            
        """
        if leica == True:
            add_leica = "-leica mdb -O.o '" + self.observer + "' -O.ag '" + \
            self.institution + "' -O.mo '" + site + "' -O.mov " + str(moving) + " "
        else:
            add_leica = ""
            # Extract information about the file.
            status = shellcmd("teqc " + add_leica + "++config " + input_file + r" | grep '\-O\.mo\[nument\]' -a")
            if status['stderr'] != '' and status['stderr'].find("Notice") == -1:
                print(status['stderr'])
                return
            ret = status['stdout']
            site = ret[15:19].lower()
            status = shellcmd("teqc " + add_leica + "++config " + input_file + r" | grep '\-O\.st\[art\]' -a")
            if status['stderr'] != '' and status['stderr'].find("Notice") == -1:
                print(status['stderr'])
                return
            start_date = status['stdout']
            
            start_date = datetime.datetime(int(start_date[11:15]),
                                           int(start_date[16:18]), 
                                           int(start_date[19:21]))
            
            status = shellcmd("teqc " + add_leica + "++config " + input_file + r" | grep '\-O\.e\[nd\]' -a")
            if status['stderr'] != '' and status['stderr'].find("Notice") == -1:
                print(status['stderr'])
                return
            end_date = status['stdout']
            end_date = datetime.datetime(int(end_date[11:15]),
                                         int(end_date[16:18]), 
                                         int(end_date[19:21]))
            
        start_doy = start_date.timetuple().tm_yday
        end_doy = end_date.timetuple().tm_yday
        
        print("Commencing windowing on " + input_file)
        print("This file begins on " + start_date.strftime("%Y.%m.%d (DOY %j)"))
        print("and ends on " + end_date.strftime("%Y.%m.%d (DOY %j)"))
        
        cal_date = start_date - datetime.timedelta(days=1)
        for doy in range(start_doy,end_doy + 1):
            print("Processing day " + str(doy))
            
            cmd = "teqc " + add_leica + "-st " + cal_date.strftime("%Y%m%d") + str(st_offset) + " +dh " + \
            str(dh) + " " + input_file + " > " + site + "_" + str(doy).zfill(3) + "0_ol." + \
            cal_date.strftime("%y") + "o"
            print("    " + cmd)
            status = shellcmd(cmd)
            if status['stderr'] != 'None' and status['stderr'].find("Notice") == False:
                print("         " + status['stderr'])
            # Increment date.
            cal_date = cal_date + datetime.timedelta(days=1)
        print("Done.")
        
        
       
class Kinematic:
    """Kinematic GPS processing functionality, utilising Track and Gamit."""
    
    def __init__(self):
        self.ion_stats = None
        self.MW_WL = None
        self.LG = None
        self.apriori = None
        self.config_subfolder = "gps_config/"
    

    def get_orbits(self, year, start_doy, end_doy, 
        clearup=True,
        rename=False):
        """Download daily IGS sp3 orbit files and overlap them.
        
        Each resulting overlapped file contains the previous day, the current
        day, and the next day.
        
        Do not try to download orbits spanning two years - only run on a 
        per-yearly basis.
        
        If clearup=True, the un-overlapped files will be deleted afterwards and
        overlapped files moved into the main directory. Otherwise, they will
        remain in cat_sp3/.

        If rename=True, files will have the year signifier removed from their
        filename.

        """
        start_doy = int(start_doy)
        end_doy = int(end_doy)
        n_days = end_doy - start_doy + 1 
        if n_days < 0:
            print("It looks like you've entered the number of days to download, not the end day. You need to specify the end day. Exiting...")
            return
        shellcmd("mkdir sp3_dl")
        cmd = "cd sp3_dl ; sh_get_orbits -orbit igsf -yr " + str(year) + " -doy " + \
            str(start_doy) + " -ndays " + str(n_days) + " -nofit"
        print(cmd)
        status = shellcmd(cmd)
        print(status['stdout'])
        print(status['stderr'])
        
        status = shellcmd("ls sp3_dl/*.sp3")
        files = status['stdout'].split("\n")
        shellcmd("mkdir cat_sp3")   
        
        doy = start_doy
        for prev,item,nex in neighborhood(files):
            if prev == None:
                prev = ""
            if nex == None:
                nex = ""
            if item == "":
                break
            if rename:
                newfn = "cat_sp3/igs" + str(doy).zfill(3) + ".sp3"
            else:
                yr = int(year) % 1000
                newfn = "cat_sp3/igs" + str(yr) + str(doy).zfill(3) + ".sp3"
            try:
                fn = open(newfn)
                fn.close()
                print("File for doy " + str(doy) + "already exists, skipping")
                continue
            except IOError:
                cmd = "cat " + prev + " " + item + " " + nex + \
                " > " + newfn
                print(cmd)
                shellcmd(cmd)  
            doy = doy + 1
        
        if clearup == True:
            print("Clearing up...")
            shellcmd("rm -r sp3_dl")
            shellcmd("mv cat_sp3/* .")
            shellcmd("rm -r cat_sp3")
        print("Done.")
             
    
    def crx2rnx(suffix):
        """Basic wrapper to CRX2RNX, which decompresses *.*d rinex files.

        CRX2RNX must be on system path.
        
        Inputs:
            suffix = file suffix, e.g. 11d for 2011 compressed files.
        Outputs:
            none.
        """
        files = shellcmd("ls *." + suffix)
        files = st.split(files['stout'],"\n")
        for fn in files:
            print(fn)
            shellcmd("CRX2RNX " + fn)
        print("Done")
         

    def track(self, base, rover, 
        year, doy_start, doy_end, 
        show_plot=True, 
        use_auto_qa=True, spearman_threshold=None, rms_threshold=20):
        """Wrapper to track kinematic processing.
        
        Processing takes one of two slightly different approaches. With 
        use_auto_qa=True, each day of data will be tested for linearity in 
        East versus North using Spearman correlation. If the result of the test
        is > spearman_threshold, the day is accepted automatically and 
        processing moves on to the next day with no interaction necessary. 
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
              exceed for the day to be approved automatically.
            rms_threshold: value in mm that median RMS must not exceed for day to
            be approved automatically.
        Outputs:
            None.
            
        """
                
        # Set up logging. 
        lfn = "gps.track." + rover + ".log"
        logging.basicConfig(filename=lfn,level=logging.INFO,format='%(message)s')
        logging.info("\n\n" + datetime.datetime.now().strftime("%Y-%m-%d %I:%M:%S %p") +
        ": BEGINNING NEW PROCESSING BATCH")
        logging.info("Start: day " + str(doy_start) + ", End: day " + \
        str(doy_end) + ", Base/Static: " + base)

        yr_short = year % 1000
        
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
                    ion_stats = input("    Ion Stats: ")
                    if len(ion_stats) == 0: ion_stats = self.ion_stats
                    MW_WL = input("    MW_WL Weighting: ")
                    if len(MW_WL) == 0: MW_WL = self.MW_WL
                    LG = input("    LG Combination Weighting: ")
                    if len(LG) == 0: LG = self.LG
                    print("Processing with new parameter values...")
                else:
                    print("Processing with defaults...")
                        
                # Construct track argument.
                # Parameters (-s) are in exact order expected by the cmd file, 
                # do not change!!
                outf = rover + "_" + base + "_" + str(year) + "_" + str(doy).zfill(3) + ".out"
                
                cmd = 'track -f {cf} -d {d} -s {ap} {mw} {lg} {istat} {base} {rover} {yr} > {out}'.format(
                    cf=os.path.join(self.config_subfolder, 'track_%s.cmd' %base),
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
                self.view_track_output(base, rover, doy, 
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
                # Ensure the directory is available.
                if os.path.isdir('processed') == False:
                    print("Making processed/ subdirectory.")
                    shellcmd("mkdir processed")
                # Move and rename NEU, GEOD results
                shellcmd("cp track.{ftype}.{r}.LC processed/{r}_{b}_{y}_{d}_{ftype}.dat".format(ftype='NEU', **save_opts))
                shellcmd("cp track.{ftype}.{r}.LC processed/{r}_{b}_{y}_{d}_{ftype}.dat".format(ftype='GEOD', **save_opts))
                # Move figure file of results
                shellcmd("mv {f} processed/".format(f=save_plot_to))
            
        print("Batch finished.")
    
    
    def view_track_output(self, base, rover, doy, gtype='NEU', fname=None, 
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
                    ms=5, c='tab:blue', marker='x', markeredgewidth=1, linestyle='none',
                    label='Position')
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


# ------------------------------------------------------------------------------
# Command line access functionality
# This should probably be re-written to use optparse...  
if __name__ == "__main__":
    
    import sys
    k = Kinematic()
    if len(sys.argv) == 1:
        print ("""
        This is the command line interface to gps. Examine the module
        code (gps.py) for full information. 
        
        Usage:
         view_track_output [base] [rover] [doy] [file_type(default=NEU)]
         get_orbits [year] [start doy] [end doy]
         crx2rnx [suffix, e.g. 11d] 
         
        """    )     
    elif sys.argv[1] == 'view_track_output':
        if len(sys.argv) == 6:
            k.view_track_output(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
        else:
            k.view_track_output(sys.argv[2],sys.argv[3],sys.argv[4])                
    elif sys.argv[1] == 'get_orbits':
        k.get_orbits(sys.argv[2],sys.argv[3],sys.argv[4])
    elif sys.argv[1] == 'crx2rnx':
        k.crx2rnx(sys.argv[2].strip()) 
    else:
        print ("Unknown command. Maybe the function hasn't been implemented for command line access?")
                
                
