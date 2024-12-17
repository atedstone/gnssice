""" 
Download and organise IGS IONEX files. 

Andrew Tedstone (andrew.tedstone@unifr.ch), July 2022
"""
from __future__ import annotations
import os
import subprocess
import argparse

#######
IONEX_PROVIDER = 'PL' # CDDIS, PL
IONEX_SOURCE = 'ftp://igs-final.man.olsztyn.pl/pub/gps_data/GPS_IONO/cmpcmb/'
#PL: 'https://cddis.nasa.gov/archive/gnss/products/ionex/'
# see https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html
#######

def download(
    year : str,
    doy_start : int,
    doy_finish : int,
    dl_location: str=IONEX_SOURCE
    ) -> None:
    """
    Download IONEX files to current working directory.
    """
    
    yr = str(year % 1000).zfill(2)

    for doy in range(doy_start, doy_finish):
        d = str(doy).zfill(3)
        if IONEX_PROVIDER == 'CDDIS':
            fpath = '{yr}/{d}/igsg{d}0.{yr}i.Z'.format(yr=yr, d=d)
        elif IONEX_PROVIDER == 'PL':
            #fpath = '{yr}{d}/igsg{d}0.{yr}i.Z'.format(yr=yr, d=d)
            fpath = '{yr}{d}'.format(yr=yr, d=d)
        else:
            raise ValueError('Unknown IONEX_PROVIDER')

        # There was a change of filename format in 2023...
        # https://files.igs.org/pub/resource/guidelines/Guidelines_For_Long_Product_Filenames_in_the_IGS_v2.0.pdf
        print(yr, doy)
        if ((int(yr) == 23 and int(doy) > 42) or (int(yr) > 23)) and IONEX_PROVIDER == 'PL':
            print('in here')
            # AAAVPPPTTT_YYYYDDDHHMM_LEN_SMP_[SSSSMRCCC_]CNT.FMT[.gz]
            # IGS0OPSFIN_YYYYDDD0000_01D_02H_GIM.INX.gz
            fpath = os.path.join(fpath,'IGS0OPSFIN_{y}{d}0000_01D_02H_GIM.INX.gz'.format(y=year, d=d))
            rename = 'igsg{d}0.{yr}i.Z'.format(yr=yr, d=d)
        else:
            fpath = os.path.join(fpath, 'igsg{d}0.{yr}i.Z'.format(yr=yr, d=d))
            rename = False

        def mkcmd(fpath, rename):
            fullp = os.path.join(dl_location, fpath)
            if rename is not False:
                cmd = 'wget -cO - %s > %s' %(fullp, rename)
            else:
                cmd = 'wget %s' %fullp
            print(cmd)
            return cmd
        subprocess.check_output(mkcmd(fpath, rename), shell=True)
    

def unzip() -> None:
    """ Unzip IONEX files in current working directory. """
    subprocess.check_output('gunzip *.Z', shell=True)


def overlap(
    yr : str, 
    doy_start : int, 
    doy_finish : int
    ) -> None:
    """ Overlap IONEX files in current working directory """

    for doy in range(doy_start+1, doy_finish):
        newfn = 'igsg{d}0_ol.{y}i'.format(d=str(doy).zfill(3), y=yr)

        if not os.path.exists(newfn):
            fmt = 'igsg{d}0.{y}i'
            prev = fmt.format(d=str(doy-1).zfill(3), y=yr)
            curr = fmt.format(d=str(doy).zfill(3), y=yr)
            nex = fmt.format(d=str(doy+1).zfill(3), y=yr)
            cmd = 'cat {p} {c} {n} > {out}'.format(p=prev, c=curr, n=nex, out=newfn)
            print(cmd)
            subprocess.check_output(cmd, shell=True)


def cli():

    p = argparse.ArgumentParser('Download and organise IGS IONEX files in current working directory.')

    p.add_argument('year', type=int)
    p.add_argument('doy_start', type=int, help='Day of year to start on')
    p.add_argument('doy_finish', type=int, help='Day of year to finish on')


    p.add_argument('-do', type=str, choices=['all', 'download', 'unzip', 'overlap'], default='all')
    args = p.parse_args()

    yr = str(args.year % 1000).zfill(2)

    if args.do == 'all' or args.do == 'download':
        print('Downloading ...')
        print('Provider:', IONEX_PROVIDER)
        download(args.year, args.doy_start, args.doy_finish)
        print('\n')

    if args.do == 'all' or args.do == 'unzip':
        print('Unzipping ...')
        unzip()
        print('\n')

    if args.do == 'all' or args.do == 'overlap':
        print('Overlapping ...')
        overlap(yr, args.doy_start, args.doy_finish)
        print('\n')



