#!/usr/bin/env python3
"""
Script to download VMF3 files for use by TRACK.

Andrew Tedstone, July 2025. Written using ChatGPT.
"""

import argparse
import datetime
import os
import requests
from pathlib import Path
from tqdm import tqdm

BASE_URL = "http://vmf.geo.tuwien.ac.at/trop_products/GRID/1x1/VMF3/VMF3_OP/"
HOURS = ["00", "06", "12", "18"]

def download_file(url, output_path):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        return True
    return False

def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days) + 1):
        yield start_date + datetime.timedelta(n)

def main():
    parser = argparse.ArgumentParser(
        description="Download VMF3 files for a given date range"
    )
    parser.add_argument(
        "--start-date", required=True, help="Start date in YYYY-MM-DD format"
    )
    parser.add_argument(
        "--end-date", required=True, help="End date in YYYY-MM-DD format"
    )
    parser.add_argument(
        "--output-dir", default="vmf", help="Directory to save downloaded files"
    )

    args = parser.parse_args()

    try:
        start_date = datetime.datetime.strptime(args.start_date, "%Y-%m-%d").date()
        end_date = datetime.datetime.strptime(args.end_date, "%Y-%m-%d").date()
    except ValueError:
        print("Error: Dates must be in YYYY-MM-DD format.")
        return

    if start_date > end_date:
        print("Error: Start date must not be after end date.")
        return

    output_dir = Path(args.output_dir)
    session = requests.Session()

    for date in tqdm(list(daterange(start_date, end_date)), desc="Dates"):
        year_str = date.strftime("%Y")
        ymd = date.strftime("%Y%m%d")
        for hh in HOURS:
            filename = f"VMF3_{ymd}.H{hh}"
            url = f"{BASE_URL}{year_str}/{filename}"
            output_path = output_dir / year_str / filename
            if output_path.exists():
                tqdm.write(f"Skipping (already exists): {output_path}")
                continue
            tqdm.write(f"Downloading: {url}")
            success = download_file(url, output_path)
            if not success:
                tqdm.write(f"Failed to download: {url}")

if __name__ == "__main__":
    main()
