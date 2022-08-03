#!/usr/bin/env python
""" 
Download and organise IGS IONEX files. 

Andrew Tedstone (andrew.tedstone@unifr.ch), July 2022
"""
from __future__ import annotations
import os
import subprocess
import argparse

#######
IONEX_SOURCE = 'ftp://igs-final.man.olsztyn.pl/pub/gps_data/GPS_IONO/cmpcmb/'
#######

def download(
    yr : str,
    doy_start : int,
    doy_finish : int,
    dl_location: str=IONEX_SOURCE
    ) -> None:
    """
    Download IONEX files to current working directory.
    """
    
    for doy in range(doy_start, doy_finish):
        d = str(doy).zfill(3)
        fpath = '{yr}{d}/igsg{d}0.{yr}i.Z'.format(yr=yr, d=d)
        fullp = os.path.join(dl_location, fpath)
        cmd = 'wget %s' %fullp
        print(cmd)
        subprocess.check_output(cmd, shell=True)


def unzip() -> None:
    """ Unzip IONEX files in current working directory. """
    subprocess.check_output('gunzip *.Z')


def overlap(
    yr : str, 
    doy_start : int, 
    doy_finish : int
    ) -> None:
    """ Overlap IONEX files in current working directory """

    for doy in range(doy_start+1, doy_finish):
        newfn = 'igsg{d}0_ol.{y}i'.format(d=str(doy).zfill(3), y=yr)

        if not os.path.exists(newfn):
            fmt = 'igsg{d}.{y}i'
            prev = fmt.format(d=str(doy-1).zfill(3), y=yr)
            curr = fmt.format(d=str(doy).zfill(3), y=yr)
            nex = fmt.format(d=str(doy+1).zfill(3), y=yr)
            cmd = 'cat {p} {c} {n} > {out}'.format(p=prev, c=curr, n=nex, out=newfn)
            print(cmd)
            subprocess.check_output(cmd, shell=True)


if __name__ == '__main__':

    p = argparse.ArgumentParser('Download and organise IGS IONEX files in current working directory.')

    p.add_argument('year', type=int)
    p.add_argument('doy_start', type=int, help='Day of year to start on')
    p.add_argument('doy_finish', type=int, help='Day of year to finish on')

    p.add_argument('-do', type=str, action='extend', nargs='*', default=['download', 'unzip', 'overlap'])

    args = p.parse_args()

    yr = str(args.year % 1000).zfill(2)

    if 'download' in args.do:
        print('Downloading ...')
        download(yr, args.doy_start, args.doy_finish)
        print('\n')

    if 'unzip' in args.do:
        print('Unzipping ...')
        unzip()
        print('\n')

    if 'overlap' in args.do:
        print('Overlapping ...')
        overlap(yr, args.doy_start, args.doy_finish)
        print('\n')



