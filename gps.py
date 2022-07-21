#!/usr/bin/python
"""Kinematic GPS Processing

Classes:
RinexConvert -- functionality for converting raw leica mdb files to the 
research group rinex file specification. Wrapper to TEQC.
Kinematic -- wrapper to kinematic GPS processing using Track.
PostProcess -- concatenation of Track output files, conversion to NEU, etc.

Functions:
shellcmd -- used to run subprocesses through shell.
neighborhood -- search behind/current/ahead in list simulataneously (generator).
confirm -- Prompt for yes/no response from user.

Pre-requisites:
As well as all required module imports, TEQC/Gamit/Track must be available
in the environment.

More information:
See the instruction file gps_python_instructions.txt.

OUTSTANDING ISSUES (at 15 Nov 2012):
    PostProcess.concatenate_GEOD -- if using appendToFile functionality, the 
    code currently WILL NOT skip past datasets listed in concatenate that are
    already in appendToFile - so make sure these are commented out in the XML 
    config file for the rover.
    
    Can be problems with producing NEU plots after concatenate_GEOD has run.

HISTORY
Created on Thu Feb 02 10:50:25 2012 
2022-04: Upgrade to Py3. Some refactoring to Pandas.

@author: Andrew Tedstone (andrew.tedstone@unifr.ch)

"""
# try:
#     import matplotlib
#     #matplotlib.use('GTkAgg')
#     import matplotlib.pyplot as plt
# except:
#     print("--------------------------------------------------------------")
#     print("ATTENTION: matplotlib could not load. No figures can be drawn.")
#     print("Did you start XMing? Does PuTTY have X11 forwarding enabled?"  )
#     print("--------------------------------------------------------------")
import matplotlib.pyplot as plt

import numpy as np
import scipy.stats
import subprocess
import datetime
import calendar as cal
import string as st
import math
import logging
import os
import xml.etree.ElementTree as etree


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
                                stderr=subprocess.PIPE, text=True)
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
    item = iterator.next()  # throws StopIteration if empty.
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
            rinex_st = shellcmd("teqc ++config " + fn + " | grep '\-O\.st\[art\]' -a")
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
                                 start_date=[2011,5,1],end_date=[2011,8,15])
        E.g. Process from a rinex file:
            window_overlap("rinexfile.11o","-22:00:00",28)
        
        
        Outputs:
            None.
            
        """
        if leica == True:
            add_leica = "-leica mdb -O.o '" + self.observer + "' -O.ag '" + \
            self.institution + "' -O.mo '" + site + "' -O.mov " + str(moving) + " "
            start_date = datetime.datetime(start_date[0], start_date[1],
                                           start_date[2])
            end_date = datetime.datetime(end_date[0], end_date[1], 
                                         end_date[2])
        else:
            add_leica = ""
            # Extract information about the file.
            status = shellcmd("teqc " + add_leica + "++config " + input_file + " | grep '\-O\.mo\[nument\]' -a")
            if status['stderr'] != '' and status['stderr'].find("Notice") == -1:
                print(status['stderr'])
                return
            ret = status['stdout']
            site = ret[15:19].lower()
            status = shellcmd("teqc " + add_leica + "++config " + input_file + " | grep '\-O\.st\[art\]' -a")
            if status['stderr'] != '' and status['stderr'].find("Notice") == -1:
                print(status['stderr'])
                return
            start_date = status['stdout']
            
            start_date = datetime.datetime(int(start_date[11:15]),
                                           int(start_date[16:18]), 
                                           int(start_date[19:21]))
            
            status = shellcmd("teqc " + add_leica + "++config " + input_file + " | grep '\-O\.e\[nd\]' -a")
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
    

    def get_orbits(self,year,start_doy,end_doy,clearup=True):
        """Download daily IGS sp3 orbit files and overlap them.
        
        Each resulting overlapped file contains the previous day, the current
        day, and the next day.
        
        Do not try to download orbits spanning two years - only run on a 
        per-yearly basis.
        
        If clearup=True, the un-overlapped files will be deleted afterwards and
        overlapped files moved into the main directory. Otherwise, they will
        remain in cat_sp3/.
        """
        start_doy = int(start_doy)
        end_doy = int(end_doy)
        n_days = end_doy - start_doy + 1 
        if n_days < 0:
            print("It looks like you've entered the number of days to download, not the end day. You need to specify the end day. Exiting...")
            return
        shellcmd("mkdir sp3_dl")
        cmd = "cd sp3_dl ; sh_get_orbits -orbit igsf -archive sopac -yr " + str(year) + " -doy " + \
            str(start_doy) + " -ndays " + str(n_days) + " -nofit"
        print(cmd)
        status = shellcmd(cmd)
        print(status['stdout'])
        print(status['stderr'])
        
        status = shellcmd("ls sp3_dl/*.sp3")
        files = st.split(status['stdout'],"\n")
        shellcmd("mkdir cat_sp3")   
        
        doy = start_doy
        for prev,item,nex in neighborhood(files):
            if prev == None:
                prev = ""
            if nex == None:
                nex = ""
            if item == "":
                break
            newfn = "cat_sp3/igs" + str(doy).zfill(3) + ".sp3"
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
        doy_start, doy_end, 
        show_plot=True, 
        use_auto_qa=True, spearman_threshold=0.66):
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
            doy_start: day of year to start on.
            doy_end: day of year to end on.
            show_plot: if True, plot of data will be popped up automatically.
            use_auto_qa: if True, high quality days will automatically be 
              accepted. When a day is automatically accepted, a plot of the 
              data will be saved but will not be displayed on screen, even 
              if show_plot=True.
           spearman_threshold: the value that the spearman coefficient must 
              exceed for the day to be approved automatically.
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
        
        # Enter main processing loop, works on a per-day basis
        for doy in range(doy_start,doy_end):
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
                status = shellcmd("convertc " + str(vals[3]) + " " + 
                str(vals[4]) + " " + str(vals[5]) + " XYZ")                
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
                    exclude_svs = input("    Exclude satellites, if multiple separate by single space (exclude_svs): ")
                    if len(exclude_svs) == 0: exclude_svs = ''
                    print("Processing with new parameter values...")
                else:
                    print("Processing with defaults...")
                        
                # Construct track argument.
                # Parameters (-s) are in exact order expected by the cmd file, 
                # do not change!!
                outf = rover + "_" + base + "_" + str(doy).zfill(3) + ".out"
                
                cmd = "track -f " + self.config_subfolder + "track_" + base + ".cmd -d " + str(doy).zfill(3) + " -s " + \
                apriori + " " + str(MW_WL) + " " + str(LG) + " " + \
                str(ion_stats) + " " + base + " " + rover + " " + str(exclude_svs) + " > " + \
                outf
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
                fid.close()
                if track_error == True:
                    break # break out of while:true loop
                
                # Do automated quality check, if requested.
                if use_auto_qa == True: 
                    data = self.read_track_file("track.NEU." + rover + ".LC")                      
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
                
                # If day passed automatically, prevent the plot from popping up
                if use_auto_qa == True and keep == True:
                    __show_plot = False
                # Otherwise continue to use user-specified option.
                else:
                    __show_plot = show_plot
                    
                # Save a scatter plot, also display subject to above.
                plot_fname = "track.NEU." + rover + ".LC"
                ret_fname = self.view_track_output(base, rover, doy, 
                                       fname=plot_fname,
                                       display=__show_plot)
                
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
                    comment = 'Spearman: ' + str(spearman[0])
                
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
                # Move and rename NEU results
                shellcmd("cp track.NEU." + rover + ".LC processed/" + rover + "_" +
                base + "_" + str(doy).zfill(3) + "NEU.dat")
                # Move and rename GEOD results
                shellcmd("cp track.GEOD." + rover + ".LC processed/" + rover + "_" +
                base + "_" + str(doy).zfill(3) + "GEOD.dat")
                # Move figure file of results
                shellcmd("mv " + ret_fname + " processed/")
            
        print("Batch finished.")
                
    
    def read_track_file(self, fname):
        """Reads the data from track files output in either NEU or GEOD.
        
        Inputs:
            fname - the path and name of the file to read in. 
        Outputs:
            Returns data - a dictionary of lists, keyed by column name.
            
        N.b. NotFK provides the letter identifier associated with NotF. In matlab
        these two 'columns' seem either to appear as one, or the letter identifier
        doesn't get read in at all - not sure which.... (AJT 8/Feb/2012)
        
        """   
        # Define dataset structure based on georeferencing type.
        if fname.count('NEU') > 0:
            keys = ['YY','MM','DD','HR','MIN','Sec','dNorth','dNorth+-','dEast',
                    'dEast+-','dHeight','dHeight+-','RMS','#DD','Atm','Atm+-',
                    'FractDOY','Epoch','#BF','NotF','NotFK','Rho_UA']
            converters = ([float] * 20) + [str] + [float]
        elif fname.count('GEOD'):
            keys = ['YY','DOY','Seconds','Latitude','Longitude','Height','SigN',
                    'SigE','SigH','RMS','#','Atm','Atm+-','Fract DOY','Epoch',
                    '#BF','NotF','NotFK']
            converters = ([float] * 17) + [str]
        else:
            raise ValueError('Georeferencing not of accepted types GEOD or NEU.')
    
        # create the data structure
        data = {}    
        for k in keys:
            data[k] = []     
        
        # Open the file
        try:  
            infile = open(fname,'r')
        except IOError as e:
            raise e
        # Read past column names as Python doesn't read them correctly
        infile.readline()
        # Read past the measurement units
        infile.readline()  
        
        # loop over the remaining lines in the file
        for line in infile.readlines():
            # Skip line if it contains bad data
            if line.find("*") == -1:
                line = line.split()
                for i in range(0,len(line)):
                    try:
                        data[keys[i]].append(converters[i](line[i]))
                    except ValueError as detail:
                        print (detail, ' filling with value -100')
                        data[keys[i]].append(-100)
        infile.close()
            
        return data
    
    
    def view_track_output(self, base, rover, doy, gtype='NEU', fname=None, 
                          display=True):
        """Display a scatter plot of reconciled daily track data.
        
        Inputs:
            base: 4-character identifier of base (static) receiver
            rover: 4-character identifier of moving receiver
            doy: day of year to examine
            gtype: georeferencing type. Optional. Currently only NEU supported
            fname: filename of file. Optional. If not used, the filename will
                be deduced from the other input parameters.
        
        Outputs:
            None.
        
        """
        if fname == None:        
            fname = rover + '_' + base + '_' + str(doy) + gtype + '.dat'
        
        if gtype=='NEU':    
            data = self.read_track_file(fname)
        else:
            print( 'Files georeferenced in a format other than NEU are currently unsupported.')
            return
        
        #plt.ioff()
        plt.plot(data['dEast'], data['dNorth'],
                    ms=5, c='b', marker='x', markeredgewidth=1, linestyle='none',
                    label='Position')
        plt.hold(True)
        plt.plot(data['dEast'][0],data['dNorth'][0],
                 c='r', marker='+', markersize=30,
                 label='Start')
        plt.plot(data['dEast'][-1],data['dNorth'][-1],
                 c='k', marker='+', markersize=30,
                 label='End')
        
        plt.ylabel('North (m)')
        plt.xlabel('East (m)')
        plt.grid(True)
        #plt.title('n = ' + `len(data['dEast'])`)
        plt.axes().set_aspect('equal','datalim')    
        plt.legend(numpoints=1, scatterpoints=1, markerscale=0.6, loc='best')    
        plt.hold(False)
        fname = 'trackpy.NEU.' + rover + '.LC.' + str(doy) + '.eps'
        plt.savefig(fname,
                    orientation='landscape', dpi=150)  
        if display == True:
           plt.show()
        return fname
              


class PostProcess:
    """GPS post-processing class."""
    
    def __init__(self):
        # Length of Earth's semi-major axis in metres (WGS84)
        self.a = 6378137.0
        # Length of Earth's semi-minor axis in metres (WGS84)
        self.b = 6356752.3142
        # First numerical eccentricity
        self.e2 = 1.0 - math.pow((self.b/self.a),2)
        self.config_subfolder = "gps_config/"
        
    
    def concatenate_GEOD(self,rover,returnData=False,plotData=False,
                         appendToFile=None):
        """Concatenate, correct and convert to NEU coordinates.
        
        !!! 12-Mar-2012: added appendToFile functionality. However, 
        this has not yet been tested.
        
        Concatenate all data for a rover into one dat file using specifications
        in the rover's config (XML) file. So, data from multiple years are joined
        together.
        
        If a concatenated file for a single year exists in the working directory
        then this is used rather than concatenating daily files.
        
        If a filename is specified in appendToFile then the function will 
        attempt to read in the data from it. Data in it will not be corrected. 
        Additional data not in the file will be concatenated and corrected 
        and a new dat file produced.
        
        Calculate north, east and up coordinates for each observation.
        
        Apply corrections listed in the config file.
        
        Save a file containing the concatenated data.
        
        Draw a graph of corrected data is plotData=True.
        
        The Fract DOY column runs continuously through multiple years,
        e.g. if there are data from 2009 and 2010 then all Fract DOYs in
        2010 will be +365.
        
        Inputs:
            rover: 4-character identifier of site to concatenate.
            returnData: if true, smap_all_save is returned.
            plotData: if true, data are plotted by plot_neu after concatenation.
            appendToFile: if a filename, function attempts read data in from it.
             Function will not apply any corrections to data already in file.
             Constituent datasets should still be listed in the config file
             as this information is used to name to output dataset.
        Outputs:
            smap_all_save: only if returnData=True.
            
            
        Useful thing to note: NEU vs XYZ is mind-bending. Here, the 'neu' 
        variable has three arrays: 'x', 'y', and 'z'. These correspond to 
        'n', 'e' and 'u' - in that order!!! So x is North. Hopefully this will
        save somebody falsely debugging the routine yet again.... AJT 27/05/13
        
        """
        # Load the XML configuration file        
        try:
            conf = etree.parse(self.config_subfolder + "config_" + rover + ".xml")
        except IOError:
            print("Configuration file does not appear to exist. Exiting.")
            return
        root = conf.getroot()
        
        # Initialise local variables        
        years_as_doy = 0  # The number of doys to add on according to years
        last_year = 0     # The year of the previous dataset processed
        start_year = 0    # First year specified for concatenation
        n_dataset = 1     # Dataset counter
        labels = []       # To pass to plotting function
        pre_res = None    # Holds the results of from file appendToFile
        last_conc_rec = None # If pre_res != None, holds year and fract doy of 
                             # last record in file appendToFile   
        pre_res_fractdoy_sorted = False 
        
        # If an already-concatenated multi-year file exists, sort it out here.       
        if appendToFile != None:
            try:
                print('Loading ' + appendToFile + '...')
                pre_res = np.loadtxt(appendToFile)
                print('...file loaded.')
            except IOError:
                print("Could not open filename specified by appendToFile.")
        # If file was opened, retrieve details of last record
        if pre_res != None:
            last_conc_rec = {"year":int(pre_res[-1,0]),"fract_doy":pre_res[-1,13]}
            # Check whether datasets already have integer identifiers.
            if pre_res.shape[1] < 21:
                ids = np.ones((pre_res.shape[0],1)) * n_dataset
                pre_res = np.concatenate((pre_res,ids),axis=1)
                n_dataset += 1
            else:
                n_dataset = pre_res[-1,20] + 1
              
        # Do concatenation
        for spec in root.find("concatenate"):
            
            # Skip trying to concatenate individual results files if a 
            # multi-year file has already been read in
#             if pre_res != None:
#                 if int(last_conc_rec["year"]) > int(spec.attrib["year"]):
#                     continue
#                 if int(last_conc_rec["fract_doy"]) > spec.attrib["from"]:
#                     if last_conc_rec["fract_doy"] < spec.attrib["to"]:
#                         print "Warning: last day in provided file falls inside the date bounds of another dataset."
#                     continue

                        
            excluded = []  # List of doys to exclude
            labels.append([n_dataset,spec.attrib['year'] + " " + spec.attrib['base'] + 
            " " + spec.attrib['extra_id']])
            
            # See if a dat file for this dataset has already been created
            try:
                if spec.attrib['extra_id'] != '':
                    eid = spec.attrib['extra_id'] + "_"
                else:
                    eid = ""
                
                # Try to load data directly into numpy array
                # Filename schema: <base>_<rover>_<year>_<extra id>geod.dat
                fname = spec.attrib['base'] + "_" + rover + "_" + \
                              str(spec.attrib['year']) + "_" + eid + \
                              "geod.dat"
                print("Seeing if " + fname + " exists."  )
                res = np.loadtxt(fname)
                
                # Check to see if it has already been processed to NEU
                if res.shape[1] == 17:
                    print ("File loaded.")
                elif res.shape[1] == 20:
                    print ("File loaded. N, E and U seem to be calculated already. They will be re-calculated...")
                    res = res[:,0:17]   
                elif res.shape[1] == 21:
                    print ("File loaded. There is a fileset ID column.  N, E and U seem to be calculated already. They will be re-calculated...")
                    res = res[:,0:17]
                else:
                    print ("File found but unexpected number of columns. Regenerating from daily files.")
                    raise(IOError)
                
            except IOError:
                print ("Reading daily files...")
           
                # Get excluded doys
                for e in spec.findall("ex"):
                    from_doy = e.attrib['from']
                    to_doy = e.attrib['to']
                    r = range(int(from_doy),int(to_doy) + 1)
                    excluded.extend(r)
                
                print ("Concatenating " + str(spec.attrib['year']) + " " + \
                spec.attrib['base'] + ", excluding days " + str(excluded))
                res = self.concatenate_daily_GEOD(spec.attrib['base'],
                                                  rover,spec.attrib['year'],
                                                  exclude_doy=excluded)
            
            # Figure out how many years worth of days to add to FractDOY
            if last_year == 0:
                # Increment from end of multiple-years file if provided
                if pre_res != None and pre_res_fractdoy_sorted == False:
                    # ensures we only do this once!
                    pre_res_fractdoy_sorted = True
                    last_year = int(last_conc_rec["year"])
                    start_year = pre_res[0,0]
                    if last_year == int(spec.attrib['year']):
                        print ('new dataset is in same year as pre_res last record.')
                        upto_year = last_year - 1
                    else:
                        upto_year = last_year
                    for y in range(int(start_year),int(upto_year+1)):
                        if cal.isleap(y):
                            years_as_doy += 366
                        else:
                            years_as_doy += 365
                else:
                    last_year = int(spec.attrib['year'])
                    # We need start_year to assist with corrections later
                    start_year = int(spec.attrib['year'])
            # This takes account of potential for multiple datasets in one year
            elif last_year < int(spec.attrib['year']):
                last_year = int(spec.attrib['year'])
                if cal.isleap(int(spec.attrib['year'])):
                    years_as_doy = years_as_doy + 366
                else:
                    years_as_doy = years_as_doy + 365
                        
            # Concatenate these results to the main results array
            if res != None:  
                # Add DOYs on    
                res[:,13] = res[:,13] + years_as_doy
                # Add a dataset identity column onto the end
                ids = np.ones((res.shape[0],1)) * n_dataset
                res = np.concatenate((res,ids),axis=1)
                n_dataset += 1
                if 'smap_all' not in locals():    
                    smap_all = res
                else:
                    smap_all = np.concatenate((smap_all,res),axis=0)
                    
        # Sort by fractional DOY    
        smap_all = smap_all[smap_all[:,13].argsort(),]
        
        # Convert elliposidal coordinates to cartesian       
        tmp_lat = smap_all[:,3] * (math.pi/180)
        tmp_lon = smap_all[:,4] * (math.pi/180)
        xyz = self.ell2xyz(tmp_lat,tmp_lon,smap_all[:,5]) # 5 is height
        
        # Convert ECEF XYZ to NEU
        if pre_res != None:
            # We have to recalculate local cartesian for pre_res first 100 points as these are not saved in pre_res.
            pre_res_xyz = self.ell2xyz(pre_res[0:100,3] * (math.pi/180),
                                       pre_res[0:100,4] * (math.pi/180),
                                       pre_res[0:100,5]) 
            dX = xyz['x'] - np.median(pre_res_xyz['x']) 
            dY = xyz['y'] - np.median(pre_res_xyz['y'])
            dZ = xyz['z'] - np.median(pre_res_xyz['z'])
            # The calculated neu coordinates reference lat-lon origin point.
            # All the nasty matrix reshaping is to satisfy the requirements of
            # the function!
            neu = self.ct2lg(np.array([dX]).T,np.array([dY]).T,np.array([dZ]).T,
                       np.array([pre_res[0,3] * (math.pi / 180)]),
                       np.array([pre_res[0,4] * (math.pi / 180)])) 
                       
        else:
            dX = xyz['x'] - np.median(xyz['x'][0:100]) 
            dY = xyz['y'] - np.median(xyz['y'][0:100])
            dZ = xyz['z'] - np.median(xyz['z'][0:100])
            # The calculated neu coordinates reference lat-lon origin point.
            # All the nasty matrix reshaping is to satisfy the requirements of
            # the function!
            neu = self.ct2lg(np.array([dX]).T,np.array([dY]).T,np.array([dZ]).T,
                           np.array([smap_all[0,3] * (math.pi / 180)]),
                           np.array([smap_all[0,4] * (math.pi / 180)]))        
        
        
        # Do corrections
        pre_n = 0
        pre_e = 0
        pre_u = 0
        pre_corr_applied = False
        for corr in root.find("correct"):
            print ("Applying corrections.")
                        
            # Need to add on doys according to year
            add_on = 0
            for y in range(int(start_year),int(corr.attrib['year'])):            
                if cal.isleap(y):
                    add_on = add_on + 366
                else:
                    add_on = add_on + 365   
            beyond_doy = float(corr.attrib['doy']) + add_on
            
            # Work out corrections to apply. If these corrections are
            # within the time/date range of the data in appendToFile, don't
            # apply them (as the file is already corrected), but do add them
            # all up so that their incremental effect can be added to the first
            # uncorrected record.
            if last_conc_rec != None and pre_corr_applied == False:
                # This correction is within time of appendToFile, so just
                # increment it to total correction variables then move on
                # to next correction
                if float(corr.attrib["doy"]) + add_on < last_conc_rec["fract_doy"]:
                    pre_n = pre_n + float(corr.attrib["n"])
                    pre_e = pre_e + float(corr.attrib["e"])
                    pre_u = pre_u + float(corr.attrib["u"])
                    continue
                # This is the first correction to be within time of 'new' data,
                # so add the total correction variables to it in preparation
                # for adding it onto the data below
                elif float(corr.attrib["doy"]) + add_on >= last_conc_rec["fract_doy"]:
                    ch_n = pre_n + float(corr.attrib['n'])
                    ch_e = pre_e + float(corr.attrib['e'])
                    ch_u = pre_u + float(corr.attrib['u'])
                    # Set flag to show that we've gone beyond the date range
                    # of the appendToFile data
                    pre_corr_applied = True
                else:
                    print ("Why are we here?")
            else:
                ch_n = float(corr.attrib['n'])
                ch_e = float(corr.attrib['e'])
                ch_u = float(corr.attrib['u'])           
            
            # Identify the data beyond the specified DOY to correct
            change = (smap_all[:,13] >= beyond_doy).nonzero()
            neu['x'][change] = neu['x'][change] + ch_n
            neu['y'][change] = neu['y'][change] + ch_e
            neu['z'][change] = neu['z'][change] + ch_u
        
        # Append NEU columns to matrix
        print ("Concatenating N,E,U to main results.")
        # Moves the identity column to the right hand side of the final matrix
        smap_all_save = np.concatenate((smap_all[:,0:17], #'up to but excluding 17'
                                       np.array([neu['x']]).T,
                                       np.array([neu['y']]).T,
                                       np.array([neu['z']]).T,
                                       np.array([smap_all[:,17]]).T), 
                                       axis=1)
                                       
        # Add multi-year results on if they are available
        # Assumes they have already been corrected
        if pre_res != None:
            smap_all_save = np.concatenate((pre_res,smap_all_save),axis=0) 
                                       
        # Write to DLM file
        # Name by first and last years
        if pre_res != None:
            first_year = str(int(pre_res[0,0]))
        else:
            first_year = root.find("concatenate")[0].attrib['year']
        last_year = root.find("concatenate")[-1].attrib['year']
        if first_year == last_year:
            yrstr = str(first_year)
            titleyrstr = yrstr
        else:
            yrstr = str(first_year) + "_" + str(last_year)
            titleyrstr = str(first_year) + "-" + str(last_year)
        
        fname = rover + "_" + yrstr + "_geod.dat"
        print ("Starting to save to " + fname)
        np.savetxt(fname,smap_all_save,fmt="%.12g")
        print ("File saved.")
       
       
        if plotData == True:
            
            plot_doy = smap_all[:,13]
            plot_id = smap_all[:,17]
            
            if pre_res != None:
                neu['x'] = np.concatenate((pre_res[-100000:-1,17],neu['x']))
                neu['y'] = np.concatenate((pre_res[-100000:-1,18],neu['y']))
                neu['z'] = np.concatenate((pre_res[-100000:-1,19],neu['z']))
                plot_doy = np.concatenate((pre_res[-100000:-1,13],plot_doy))
                plot_id = np.concatenate((np.array([0]*99999),plot_id))
            
            neu['doy'] = plot_doy
            neu['id'] = plot_id
#            smap_plot = np.concatenate((np.array([smap_all[:,13]]).T,
#                                    np.array([neu['x']]).T,
#                                    np.array([neu['y']]).T,
#                                    np.array([neu['z']]).T,
#                                    np.array([smap_all[:,17]]).T),axis=1)
            titlestr = rover + " " + titleyrstr
            savestr = rover + "_" + str(first_year) + "-" + str(last_year)
            labels.append([0,'pre_res'])
            self.plot_neu(neu,labels,titlestr,savestr)
        
        if returnData == True:
            return smap_all_save

        
    
    def concatenate_daily_GEOD(self,base,rover,year,start_doy=1,end_doy=366,
                               exclude_doy=[],save_to_file=False):
        """Return one large matrix of all daily files combined within 
           specified DOY range.
           
           Inputs:
               base:        4-letter identifier of base station.
               rover:       4-letter identifier of rover (on-ice GPS stake)
               year:        Year in which observations made. 
               start_doy:   Day of year on which to start concatenation
               end_doy:     Day of year on which to end concatenation
               exclude_doy: List of days to exclude from concatenation, usually
                            because of poor data quality. Optional.
               save_to_file:If true a file will be saved with standard filename
                            format.
           Outputs:
               smap_all :   Matrix of all concatenated data. Columns are in same
                            order as GEOD daily files.
                            
        """
        k = Kinematic()
        print ("Key: E=Excluded, S=Skipped (no file), C=Concatenated.")
        for doy in range(start_doy,end_doy):
            if doy in exclude_doy:
                print ('E'+ str(doy) + ', ',)
                continue
            
            # Try to read in the daily GEOD file 
            try:      
                data = k.read_track_file(rover + '_' + base + '_' + str(doy).zfill(3) + 'GEOD.dat')
            except IOError:
                print ('S' + str(doy) + ', ',)
                continue
            print ('C' + str(doy) + ', ')        
            smap = np.array((data['YY'],data['DOY'],data['Seconds'],data['Latitude'],
                             data['Longitude'],data['Height'],data['SigN'],
                             data['SigE'],data['SigH'],data['RMS'],data['#'],
                             data['Atm'],data['Atm+-'],data['Fract DOY'],
                             data['Epoch'],data['#BF'],data['NotF']))
            
            # Remove overlaps i.e. 22:00 doy-1 --> 02:00 doy+1
            smap = smap.transpose()  
            get_rid = (smap[:,13] < float(doy)) | (smap[:,13] > float(doy) + 1) 
            smap = smap[~get_rid,:]
            
            # concatenate to smap_all        
            if 'smap_all' not in locals():
                smap_all = smap
            else:
                smap_all = np.concatenate((smap_all,smap),axis=0) 
        
        # Just in case no data was extracted        
        if 'smap_all' not in locals():
                smap_all = None
                
        if smap_all != None and save_to_file == True:
            print ("Saving to file.")
            fname = base + "_" + rover + "_" + str(year) + "_geod.dat"
            np.savetxt(fname,smap_all,fmt="%.12g")
        elif smap_all == None and save_to_file == True:
            print ("No data found. No file will be saved.")
            
        print("Done")
        return smap_all
    
    
    def plot_neu(self,dneui,labels,titlestr,savestr,display=True):
        """Plot north, east and up coordinates in separate figures.
        
        4 figures are plotted: north v. east, north v. time, east v. time,
        up v. time.
        
        Each dataset is plotted in a different colour. The labels list must
        have the same number of elements as the number of datasets supplied.
        
        Inputs:
            neui: numpy array. Columns: doy, north, east, up, dataset id.
            labels: the label/description associated with each dataset id.
            titlestr : string to add to title, e.g. lev1 2009-2011.
            savestr : filename prefix for each graph. The plot type, e.g.
            north-east, is appended by this function, as is the file type.
            display: if true, figure windows will be drawn on screen as well
            as saved. 
        Outputs:
            Returns null.
            Saves 4 figures.
        
        """ 
        print ("Plotting North, East, Up, Time.")
        print ("Unless running in IPython, each figure will have to be closed before the next will open.")
        
        # Do some sorting (makes setting axis limits easier)
        n_sort = np.sort(dneui['x']) # sorted north
        e_sort = np.sort(dneui['y']) # sorted east
        u_sort = np.sort(dneui['z']) # sorted up
        # Calculate axis limits
#        n_min = np.median(n_sort[0:100])
#        n_max = np.median(n_sort[-100:])  
#        e_min = np.median(e_sort[0:100])
#        e_max = np.median(e_sort[-100:])
#        u_min = np.median(u_sort[0:100])
#        u_max = np.median(u_sort[-100:])
        
        n = dneui['x']
        e = dneui['y']
        u = dneui['z']
        d = dneui['doy']
        i = dneui['id']
        
        # North versus East
        plt.figure(1)
        for l in labels:
            li = l[0]
            ln = l[1]
            dsr = (i == int(li)).nonzero()
            plt.plot(e[dsr],n[dsr],
                     label=ln,marker="x",ms=1,lw=0)
            plt.hold(True)
        plt.title("North versus East " + titlestr)
        plt.xlabel("East (m)")
        plt.ylabel("North (m)")
        #plt.axis([e_min,e_max,n_min,n_max])
        plt.legend(numpoints=1,scatterpoints=1,loc="best")
        plt.savefig(savestr + "North-East.png",orientation="landscape")
        plt.hold(False)
        if display == True:
            plt.show()
        
        # North versus Time
        plt.figure(2)
        for l in labels:
            li = l[0]
            ln = l[1]
            dsr = (i == li).nonzero()
            plt.plot(n[dsr],d[dsr],
                     label=ln,marker="x",ms=1,lw=0)
            plt.hold(True)
        plt.title("North versus Time " + titlestr)
        plt.ylabel("Day")
        plt.xlabel("North (m)")
        #plt.ylim(n_min,n_max)
        plt.legend(numpoints=1,scatterpoints=1,loc="best")
        plt.savefig(savestr + "North-Time.png",orientation="landscape",dpi=150)
        plt.hold(False)
        if display == True:
            plt.show()
        
        # East versus Time
        plt.figure(3)
        for l in labels:
            li = l[0]
            ln = l[1]
            dsr = (i == li).nonzero()
            plt.plot(e[dsr],d[dsr],
                     label=ln,marker="x",ms=1,lw=0)
            plt.hold(True)
        plt.title("East versus Time " + titlestr)
        plt.xlabel("East (m)")
        plt.ylabel("Day")
        #plt.xlim(e_min,e_max)
        plt.legend(numpoints=1,scatterpoints=1,loc="best")
        plt.savefig(savestr + "East-Time.png",orientation="landscape",dpi=150)
        plt.hold(False)
        if display == True:
            plt.show()
        
        # Up versus Time
        plt.figure(4)
        for l in labels:
            li = l[0]
            ln = l[1]
            dsr = (i == li).nonzero()
            plt.plot(d[dsr],u[dsr],
                     label=ln,marker="x",ms=1,lw=0)
            plt.hold(True)
        plt.title("Up versus Time " + titlestr)
        plt.xlabel("Day")
        plt.ylabel("Up (m)")
        #plt.ylim(u_min,u_max)
        plt.legend(numpoints=1,scatterpoints=1,loc="best")
        plt.savefig(savestr + "Up-Time.png",orientation="landscape",dpi=150)
        plt.hold(False)
        if display == True:
            plt.show()
            
            
       
    
    def ell2xyz(self,lat,lon,h,a=None,e2=None):
        """Convert lat lon height to local north-east-up.
        
        lat, lon and h must be numpy 1-d arrays.
        
        If a and e2 are None then their values will be obtained from the 
        class variables self.a and self.e2.
        
        Based on the matlab exchange function...:
        ELL2XYZ  Converts ellipsoidal coordinates to cartesian.
        Vectorized.
         Version: 2011-02-19
         Useage:  [x,y,z]=ell2xyz(lat,lon,h,a,e2)
                  [x,y,z]=ell2xyz(lat,lon,h)
         Input:   lat - vector of ellipsoidal latitudes (radians)
                  lon - vector of ellipsoidal E longitudes (radians)
                  h   - vector of ellipsoidal heights (m)
                  a   - ref. ellipsoid major semi-axis (m); default GRS80
                  e2  - ref. ellipsoid eccentricity squared; default GRS80
         Output:  x \
                  y  > vectors of cartesian coordinates in CT system (m)
                  z /
        
         Copyright (c) 2011, Michael R. Craymer Email: mike@craymer.com
                    
        """
        if a == None:
            a = self.a
        if e2 == None:
            e2 = self.e2
            
        v = a / np.sqrt(1 - e2 * np.sin(lat) * np.sin(lat))
        x = (v + h) * np.cos(lat) * np.cos(lon)
        y = (v + h) * np.cos(lat) * np.sin(lon)
        z = (v * (1 - e2) + h) * np.sin(lat)   
                
        toret = {}
        toret['x'] = x
        toret['y'] = y
        toret['z'] = z
        return toret
        
  
    def ct2lg(self,dX,dY,dZ,lat,lon):
        """Converts CT coordinate differences to local geodetic.
        
        All inputs must be numpy arrays.
        
        Local origin at lat,lon,h. If lat,lon are vectors, dx,dy,dz
        are referenced to origin at lat,lon of same index. If
        astronomic lat,lon input, output is in local astronomic
        system. Vectorized in both dx,dy,dz and lat,lon. See also
        LG2CT.
        Version: 2011-02-19
        Useage:  [dx,dy,dz]=ct2lg(dX,dY,dZ,lat,lon)
        Input:   dX  - vector of X coordinate differences in CT
         dY  - vector of Y coordinate differences in CT
         dZ  - vector of Z coordinate differences in CT
         lat - lat(s) of local system origin (rad); may be vector
         lon - lon(s) of local system origin (rad); may be vector
        Output:  dx  - vector of x coordinates in local system (north)
         dy  - vector of y coordinates in local system (east)
         dz  - vector of z coordinates in local system (ht)
                
        Ported to Python from original Matlab exchange version...
        Copyright (c) 2011, Michael R. Craymer 
        All rights reserved.
        Email: mike@craymer.com 

        """
        n = dX.shape[0]
        if lat.shape[0] == 1:
          lat = np.ones((n,1)) * lat
          lon = np.ones((n,1)) * lon
        R = np.zeros((3,3,n))
        
        R[0,0,:] = -np.sin(lat.T) * np.cos(lon.T)
        R[0,1,:] = -np.sin(lat.T) * np.sin(lon.T)
        R[0,2,:] = np.cos(lat.T)
        
        R[1,0,:] = -np.sin(lon.T)
        R[1,1,:] = np.cos(lon.T)
        R[1,2,:] = np.zeros((1,n))
        
        R[2,0,:] = np.cos(lat.T) * np.cos(lon.T)
        R[2,1,:] = np.cos(lat.T) * np.sin(lon.T)
        R[2,2,:] = np.sin(lat.T)
        
        RR = np.reshape(R[0,:,:],(3,n))
        dx = np.sum(RR.T * np.concatenate((dX,dY,dZ),axis=1),axis=1)
        RR = np.reshape(R[1,:,:],(3,n))
        dy = np.sum(RR.T * np.concatenate((dX,dY,dZ),axis=1),axis=1)
        RR = np.reshape(R[2,:,:],(3,n))
        dz = np.sum(RR.T * np.concatenate((dX,dY,dZ),axis=1),axis=1)
        
        toret = {}
        toret['x'] = dx
        toret['y'] = dy
        toret['z'] = dz
        return toret
        
        
    
# Command line access functionality
# This should probably be re-written to use optparse...  
if __name__ == "__main__":
#    import argparse as ap  
#    
#    parser = ap.ArgumentParser(description="Command line interface to University of Edinburgh Kinematic GPS Processing package")
#    
#    parser.add_argument("-function",metavar="-f",type=str,required=True,
#                        help="The function to run. Available: view_track_output, get_orbits, crx2rnx, concatenate_geod.")
#                        
#    parser.add_argument("-d",metavar"-doy",required=False,
#                        help="[view_track_outputs] Day of year to plot.") 
#                        
#    parser.add_argument("-base",metavar"-b",required=False,
#                        help="[view_track_outputs] Identifier of GPS base station.")
#                        
#    parser.add_argument("-rover",metavar"-r",required=False,
#                        help="[view_track_outputs, concatenate_geod] Identifier of GPS rover station.")
#                        
#    parser.add_argument("-append_from_file",metavar"-aff",required=False,
#                        help=" [concatenate_geod] The path and filename of the GEOD file to append from.")
#                        
#    parser.add_argument("-start_doy",metavar"-s",required=False,
#                        help="[get_orbits] Day of year to begin orbit download.")
#                        
#    parser.add_argument("-end_doy",metavar"-e",required=False,
#                        help="[get_orbits] Day of year to end orbit download.")
#                        
#    parser.add_argument("-file_type",metavar"-ft",required=False,
#                        help=" [view_track_output] File type for view_track_output, default is NEU.")
#                        
#    parser.add_argument("-suffix",metavar"-su",required=False,
#                        help="[crx2rnx] Suffix of file, e.g. 11d, to convert to rnx.")
    
    
    
    import sys
    k = Kinematic()
    pp = PostProcess()
    if len(sys.argv) == 1:
        print ("""
        This is the command line interface to gps. Examine the module
        code (gps.py) for full information. 
        
        Usage:
         view_track_output [base] [rover] [doy] [file_type(default=NEU)]
         get_orbits [year] [start doy] [end doy]
         crx2rnx [suffix, e.g. 11d] 
         concatenate_geod [rover] [append_from_file]
         concatenate_daily_geod [base] [rover] [year] [startdoy] [enddoy] 

         
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
    elif sys.argv[1].lower() == 'concatenate_geod':
        if len(sys.argv) == 4:
            pp.concatenate_GEOD(sys.argv[2],plotData=True,appendToFile=sys.argv[3])
        else:
            pp.concatenate_GEOD(sys.argv[2],plotData=True)
    elif sys.argv[1].lower() == 'concatenate_daily_geod':
        pp.concatenate_daily_GEOD(sys.argv[2],sys.argv[3],int(sys.argv[4]),start_doy=int(sys.argv[5]),end_doy=int(sys.argv[6]),save_to_file=True)
    else:
        print ("Unknown command. Maybe the function hasn't been implemented for command line access?")
                
                
