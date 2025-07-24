#!/usr/bin/env python3
"""
Create IONEX files which span day boundaries, for use in TRACK.

Requires env variables for IONEX locations to be available.

Andrew Tedstone, July 2025. Partially written using ChatGPT.
"""

import argparse
from datetime import datetime, timedelta
from pathlib import Path
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine daily IONEX files into a file spanning the last X hours from previous day, requested day, and first X hours of next day."
    )
    parser.add_argument("date", help="Date to process in YYYY-MM-DD format")
    parser.add_argument("hours", type=int, help="Number of hours before and after to include")
    parser.add_argument('-overwrite', action='store_true')
    return parser.parse_args()

def get_julian_day(dt):
    return dt.timetuple().tm_yday

def get_ionex_filename(dt):
    ddd = get_julian_day(dt)
    yy = dt.strftime('%y')
    f = f"igsg{ddd:03d}0.{yy}i"
    pth = os.path.join(os.environ['GNSS_PATH_IONEX_DAILY'], f)
    return pth

def read_ionex_file(path):
    with open(path, 'r') as f:
        return f.readlines()

def split_header_body(lines):
    header, body = [], []
    header_done = False
    for line in lines:
        if not header_done:
            header.append(line)
            if 'END OF HEADER' in line:
                header_done = True
        else:
            body.append(line)
    return header, body

def parse_maps(lines):
    maps = []
    current_map = []
    current_time = None
    kind = None

    for line in lines:
        if 'START OF TEC MAP' in line:
            current_map = [line]
            kind = 'TEC'
        elif 'START OF RMS MAP' in line:
            current_map = [line]
            kind = 'RMS'
        elif 'EPOCH OF CURRENT MAP' in line:
            current_map.append(line)
            y, m, d, H, M, S = map(int, line[:36].split())
            current_time = datetime(y, m, d, H, M, S)
        elif 'END OF TEC MAP' in line or 'END OF RMS MAP' in line:
            current_map.append(line)
            maps.append((current_time, kind, current_map))
            current_map = []
            current_time = None
            kind = None
        elif current_map:
            current_map.append(line)
    return maps

def filter_maps(maps, start_time, end_time):
    return [m for m in maps if start_time <= m[0] <= end_time]

def update_header(header_lines, run_date, first_epoch, last_epoch, ddd_range, year, num_tec_maps):
    new_lines = []
    for line in header_lines:
        if 'PGM / RUN BY / DATE' in line:
            prog = 'combine_ionex'.ljust(20)
            run_by = 'combine_ionex'.ljust(20)
            date_str = run_date.strftime('%-d-%b-%y %H:%M').lower().rjust(20)
            new_line = f"{prog}{run_by}{date_str}PGM / RUN BY / DATE\n"
        elif 'global ionosphere maps for day' in line.lower():
            new_line = (
                f"  global ionosphere maps for days {ddd_range[0]:03d}-{ddd_range[1]:03d}, {year:4d}                         \n"
            )
        elif 'EPOCH OF FIRST MAP' in line:
            date_str = (
                f"{first_epoch.year:5d}{first_epoch.month:3d}{first_epoch.day:3d}"
                f"{first_epoch.hour:3d}{first_epoch.minute:3d}{first_epoch.second:3d}"
            )
            new_line = f" {date_str}{' ' * 39}EPOCH OF FIRST MAP\n"
        elif 'EPOCH OF LAST MAP' in line:
            date_str = (
                f"{last_epoch.year:5d}{last_epoch.month:3d}{last_epoch.day:3d}"
                f"{last_epoch.hour:3d}{last_epoch.minute:3d}{last_epoch.second:3d}"
            )
            new_line = f" {date_str}{' ' * 39}EPOCH OF LAST MAP\n"
        elif '# OF MAPS IN FILE' in line:
            new_line = f"{num_tec_maps:6d}{' ' * 54}# OF MAPS IN FILE\n"
        else:
            new_line = line
        new_lines.append(new_line)
    return new_lines

def write_ionex_file(path, header_lines, tec_maps, rms_maps):
    with open(path, 'w') as f:
        f.writelines(header_lines)

        tec_counter = 1
        for _, _, map_lines in tec_maps:
            for line in map_lines:
                if 'START OF TEC MAP' in line:
                    new_line = f"{tec_counter:6d}{' ' * 54}START OF TEC MAP\n"
                elif 'END OF TEC MAP' in line:
                    new_line = f"{tec_counter:6d}{' ' * 54}END OF TEC MAP\n"
                else:
                    new_line = line
                f.write(new_line)
            tec_counter += 1

        rms_counter = 1
        for _, _, map_lines in rms_maps:
            for line in map_lines:
                if 'START OF RMS MAP' in line:
                    new_line = f"{rms_counter:6d}{' ' * 54}START OF RMS MAP\n"
                elif 'END OF RMS MAP' in line:
                    new_line = f"{rms_counter:6d}{' ' * 54}END OF RMS MAP\n"
                else:
                    new_line = line
                f.write(new_line)
            rms_counter += 1

def main():
    args = parse_args()
    target_date = datetime.strptime(args.date, "%Y-%m-%d")
    hours = args.hours

    out_filename = f"igsg{get_julian_day(target_date):03d}0.{target_date.strftime('%y')}i"
    out_pth = os.path.join(os.environ['GNSS_PATH_IONEX_OVERLAP'], out_filename)
    if os.path.exists(out_pth) and not args.overwrite:
        print(args.date + ': overlapped file already exists')
        return

    prev_date = target_date - timedelta(days=1)
    next_date = target_date + timedelta(days=1)

    run_date = datetime.now()

    files = {
        'prev': Path(get_ionex_filename(prev_date)),
        'curr': Path(get_ionex_filename(target_date)),
        'next': Path(get_ionex_filename(next_date)),
    }

    for label, path in files.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing required IONEX file: {path}")

    all_maps = []
    for dt, label in zip([prev_date, target_date, next_date], ['prev', 'curr', 'next']):
        lines = read_ionex_file(files[label])
        header, body = split_header_body(lines)
        maps = parse_maps(body)
        all_maps.extend(maps)
        if label == 'curr':
            curr_header = header

    start_time = target_date - timedelta(hours=hours)
    end_time = target_date + timedelta(hours=24+hours)

    filtered_maps = filter_maps(all_maps, start_time, end_time)

    tec_maps = sorted([m for m in filtered_maps if m[1] == 'TEC'], key=lambda x: x[0])
    rms_maps = sorted([m for m in filtered_maps if m[1] == 'RMS'], key=lambda x: x[0])

    first_epoch = min(m[0] for m in filtered_maps)
    last_epoch = max(m[0] for m in filtered_maps)

    ddd_range = (get_julian_day(start_time), get_julian_day(end_time))
    year = target_date.year

    num_tec_maps = len(tec_maps)

    updated_header = update_header(
        curr_header, run_date, first_epoch, last_epoch, ddd_range, year, num_tec_maps
    )

    Path(os.environ['GNSS_PATH_IONEX_OVERLAP']).mkdir(exist_ok=True)
    write_ionex_file(out_pth, updated_header, tec_maps, rms_maps)

    print(f"Output written to: {out_filename}")

if __name__ == "__main__":
    main()
