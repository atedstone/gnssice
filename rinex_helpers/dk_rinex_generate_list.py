"""
Generate GNSS RINEX3 file URLs for download from the Danish dataforsyningen FTP server.

URL format:
  ftps://ftp.dataforsyningen.dk/GNSS/RINEX3/GRL/{year}/{doy}/{station}{country}_R_{year}{doy}0000_01D_01S_MO.crx.gz

Usage examples:
  # Single date
  python generate_gnss_filenames.py --start 2025-02-09

  # Date range
  python generate_gnss_filenames.py --start 2025-01-01 --end 2025-03-31

  # Custom station/country
  python generate_gnss_filenames.py --start 2025-02-09 --station KLSQ --country GRL

  # Save URLs to a file
  python generate_gnss_filenames.py --start 2025-01-01 --end 2025-12-31 --output urls.txt
"""

import argparse
from datetime import date, timedelta


BASE_URL = "ftps://ftp.dataforsyningen.dk/GNSS/RINEX3/{country}/{year}/{doy:03d}"
FILENAME  = "{station}00{country}_R_{year}{doy:03d}0000_01D_01S_MO.crx.gz"


def build_url(target_date: date, station: str = "KLSQ", country: str = "GRL") -> str:
    """Return the full FTP URL for a given date, station and country."""
    year = target_date.year
    doy  = target_date.timetuple().tm_yday          # day-of-year, 1-based

    directory = BASE_URL.format(country=country, year=year, doy=doy)
    filename  = FILENAME.format(station=station.upper(), country=country.upper(),
                                year=year, doy=doy)
    return f"{directory}/{filename}"


def generate_urls(start: date, end: date,
                  station: str = "KLSQ", country: str = "GRL") -> list[str]:
    """Generate one URL per day for [start, end] inclusive."""
    urls = []
    current = start
    while current <= end:
        urls.append(build_url(current, station, country))
        current += timedelta(days=1)
    return urls


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate GNSS RINEX3 FTP download URLs by date range."
    )
    parser.add_argument("--start",   required=True,
                        help="Start date in YYYY-MM-DD format")
    parser.add_argument("--end",     default=None,
                        help="End date in YYYY-MM-DD format (default: same as --start)")
    parser.add_argument("--station", default="KLSQ",
                        help="4-character station code (default: KLSQ)")
    parser.add_argument("--country", default="GRL",
                        help="3-character country/region code (default: GRL)")
    parser.add_argument("--output",  default=None,
                        help="Write URLs to this file instead of stdout")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    start = date.fromisoformat(args.start)
    end   = date.fromisoformat(args.end) if args.end else start

    if end < start:
        raise SystemExit("Error: --end must not be before --start.")

    urls = generate_urls(start, end, station=args.station, country=args.country)

    if args.output:
        with open(args.output, "w") as f:
            f.write("\n".join(urls) + "\n")
        print(f"Wrote {len(urls)} URL(s) to {args.output}")
    else:
        for url in urls:
            print(url)


if __name__ == "__main__":
    main()
