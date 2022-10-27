#!/usr/bin/env python
import datetime as dt
import argparse
import os, sys
import re

parser = argparse.ArgumentParser(description='reads filenames from command line. links file to fnnn.nc where nnn is the 3-digit forecast lead time. or if -m option, just output fcst time in minutes')
parser.add_argument("-m", '--minutes', action='store_true', help='output forecast time in minutes (no linking)')
parser.add_argument('filenames', type=str, nargs='+', help='filenames')
args = parser.parse_args()

current_directory = os.getcwd()
filenames = args.filenames

# Find /yyyymmddhh somewhere in the current working directory string
match = re.search(r"/([12]\d{3}[01]\d[0123]\d[012]\d)",current_directory)
if match:
    # Make a datetime object from the string.
    idate = dt.datetime.strptime(match.groups()[0], '%Y%m%d%H')
else:
    print("could not parse date from",current_directory)
    sys.exit(1)

for filename in filenames:
    # Find .yyyy-mm-dd_hh.mm.00 somewhere in the filename.
    match = re.search(r'\.(\d{4}-\d{2}-\d{2}_\d{2}\.\d{2}\.00)',filename)
    # Make a datetime object from the string.
    vdate = dt.datetime.strptime(match.groups()[0], '%Y-%m-%d_%H.%M.%S')
    # Sanity checks
    if vdate.minute != 0 or vdate.second != 0:
        print("non-zero minute or second")
    if vdate < idate:
        print("valid date less than init date", vdate, idate)
        sys.exit(2)
    # If -m option just output one forecast time in minutes.
    if args.minutes:
        print(int((vdate-idate).total_seconds()/60))
        sys.exit(0)

    forecast_hour = (vdate-idate).total_seconds()/3600
    link = 'f'+ "{0:03.0f}".format(forecast_hour) + ".nc"
    if os.path.exists(link):
        os.remove(link)
    os.symlink(filename, link)

