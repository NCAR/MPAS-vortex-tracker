import argparse
import atcf
import cartopy.feature
import datetime
import ibtracs
import logging
import matplotlib
matplotlib.use('Agg') # 'Agg' is non-interactive. No plt.show() command. TkAgg supports plt.show(). but can't get it to work as user mpasrt.
import matplotlib.pyplot as plt
import mpas
import numpy as np
import pandas as pd
import pdb
import re, os
import subprocess
import sys

"""
Read ATCF file
Optionally replace vitals with raw values from original mesh
Plot tracks meeting length and strength criteria.
Plots for different basins. 
"""
logger = logging.getLogger()
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='plot tropical cyclone track/intensity from ATCF file',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('atcf_file', type=argparse.FileType('r'), help='path to ATCF file to plot and possibily modify')
parser.add_argument('--basins', nargs='*', default=[], type=str, help='basin(s) to plot')
parser.add_argument('-d','--debug', action="store_true", help='print debug messages')
parser.add_argument('--diagdir', type=str, help='path to directory with original model diagnostic files. In case of ensemble, path to directory above member subdirectories i.e. do not include member subdirectory.')
parser.add_argument('--initfile', type=str, help='path to init.nc file with latCell,lonCell,nEdgesOnCell,cellsOnCell')
parser.add_argument('--initial_time', type=lambda d: datetime.datetime.strptime(d, '%Y%m%d%H'), help='initialization date/hour yyyymmddhh')
parser.add_argument('--models', nargs='*', type=str, help='only plot this model(s)')
parser.add_argument('--clobber', action="store_true", help='overwrite old file')
parser.add_argument('--origmesh', dest='origmesh', action="store_true", help='use original mesh')
parser.add_argument('--no-origmesh', dest='origmesh', action="store_false", help='do not use original mesh')
parser.set_defaults(origmesh=True)
parser.add_argument('--project', type=str, default="precip2020", help='project name, like precip2020 or hur15. Used to name image on web server.')
parser.add_argument('--wind_radii_method', type=atcf.iswind_radii_method, default='max', help='method to get wind extent from original mesh')
parser.add_argument('--min_warmcore_percent', type=float, default=25, help='minimum %% of warm core track times') # double percent sign to print single percent sign
parser.add_argument('--to_server', action="store_true", help='rsync to MMM server')
args = parser.parse_args()

debug = args.debug
if debug:
    logger.setLevel(logging.DEBUG)

logging.debug(args)

plotbasins   = args.basins
clobber      = args.clobber
diagdir      = args.diagdir
initfile     = args.initfile
min_warmcore_percent = args.min_warmcore_percent
models       = args.models
origmesh     = args.origmesh
project      = args.project 
wind_radii_method= args.wind_radii_method

# atcf may be a relative path, so use os.path.realpath() to get absolute path, which has initial date /yyyymmddhh/
atcf_file = os.path.realpath(args.atcf_file.name)
# origmeshFile is name of output atcf file that will have original mesh values. 
# Default name is original atcf file plus _origmesh suffix.
# Or in the case of Alex Kowaleski's WRF runs, _origmesh02 or _origmeshd03, depending on the WRF nest
if atcf_file.endswith('_origmesh'):
    origmesh = False
else: 
    origmeshFile = atcf_file + '_origmesh'

# Find yyyymmddhh in file path.
m = re.search("/[12]\d{9}/", atcf_file)
if m:
    init_yyyymmddhh = atcf_file[m.start()+1:m.end()-1]


# Skip if output file exists
# If ofile exists and you didn't request clobber, then stop.
if os.path.isfile(atcf_file+".png") and not clobber:
    logging.warning(f"{atcf_file}.png already exists. Use --clobber to override.")
    sys.exit(1)

df = atcf.read(atcf_file)

# A track has a unique combination of basin, cy, initial_time, and model
unique_track = ['basin', 'cy', 'initial_time', 'model']

# Skip lines not from a specifically-requested initialization time.
if args.initial_time:
    df = df[df.initial_time == args.initial_time]
    logging.debug(f"kept {len(df)} {args.initial_time} initialization times")

# Skip lines not from a specifically-requested model.
if models:
    df = df[df['model'].isin(models)]
    logging.debug(f"kept {len(df)} {models} lines.")


# If you need original mesh values, initialize diagdir and initfile.
if origmesh:
    if not diagdir:
        # Derive diagdir from atcf_file.
        m = re.search("/[12]\d{9}/", atcf_file)
        if m is None:
            logging.error(f"No diagdir supplied, and failed to derive from atcf_file: {atcf_file}")
            sys.exit(1)
        # path to diagnostic file original mesh
        diagdir = atcf_file[:m.end()]

    if not initfile:
        # Derive initfile from diagdir.
        # path to init.nc file with latCell,lonCell,etc.
        initfile = os.path.join(diagdir,"init.nc")

fineprint=""
warmcore_desc = f"warm core >= {min_warmcore_percent}% of the time"

# If one of unique_tracks columns is constant, put in title. 
title = ''
for col in unique_track:
    if df[col].value_counts().count() == 1:
        # Add this column value to title. It doesn't vary.
        title += " " + col + "= " + str(df[col].iloc[0])

# Make sure vmax > 30 knots for at least one time.
df = df.groupby(unique_track).filter(lambda x: x["vmax"].max() >=  30)
logging.info(f'{df.groupby(unique_track).ngroups} after vmax filter')

# overwater test not written. this is a placeholder 
#df = df.groupby(unique_track).filter(lambda x: True)
#print(df.groupby(unique_track).ngroups)

# If warmcore column exists make sure at least one time is warmcore or unknown.
df = df.groupby(unique_track).filter(atcf.iswarmcore)
logging.info(f'{df.groupby(unique_track).ngroups} after warmcore filter')

# Remove short tracks
df = df.groupby(unique_track).filter(lambda x: x["fhr"].nunique() > 2)
logging.info(f'{df.groupby(unique_track).ngroups} after short track filter')

# Remove long tracks (for debugging)
#df = df.groupby(unique_track).filter(lambda x: x["fhr"].nunique() < 9 )
#print(df.groupby(unique_track).ngroups, 'after long track filter')

def get_origmesh(track, initfile):
    basin, cy, initial_time, model = track.name
    # Do you need to find the original mesh values? Not if track has defined originalmeshfile . 
    if 'originalmeshfile' in track.columns and track["originalmeshfile"].all():
        logging.warning(f"{track_id} already has original mesh values.")
        return track
    if model == 'MPAS':
        if 'origmesh' in atcf_file:
            logging.warning(f'get_origmesh(): {track.name} already origmesh because "origmesh" is in atcf_file.')
            # TODO: assert 'origmeshfile' in track.columns or in userdefine1-4 
            return track
        logging.info(f'get raw mesh vmax and minp for {len(track.lon)} times, track {track.name}')
        # Return initfile too if you want to save the lat/lon information as a dictionary and speed things up.
        # Avoid opening and re-opening the same file over and over again.
        og, initfile = mpas.origmesh(track, initfile, diagdir, wind_radii_method=wind_radii_method, debug=debug)
    elif (model == 'ECMF' or re.compile("EE[0-9][0-9]$").match(model)):
        if all(track.rad.astype('int') == 0):
            logging.info(f"{model} getting radius of max wind and radii of 34/50/64 kt wind from {diagdir}")
            og = atcf.origgrid(track, diagdir, ensemble_prefix=["PF","CF","EE","EE0"], wind_radii_method=wind_radii_method, debug=debug)
        else:
            logging.error("That is odd. You requested original mesh values but the atcf file already has non-zero radii")
            sys.exit(1)
    elif "kowaleski/" in diagdir:
        # Tack grid nest number to end of original mesh output filename. 
        # Might help prevent applying d03 grid to d02 tracks, or vice versa.
        grid="d03"
        origmeshFile += grid
        og = atcf.origgridWRF(track, diagdir, grid=grid, wind_radii_method=wind_radii_method, debug=debug)
    else:
        logging.warning(f"no origmesh for {model}")
        return track # used to exit return code=3 but why do that if you have OFCL forecast in same file as model forecasts? that's fine.
    return og # replace track with og DataFrame so when track is plotted below, it is the original mesh stuff.

if origmesh:
    if os.path.isfile(origmeshFile) and not clobber:
        logging.warning(f"{origmeshFile} already exists. Use --clobber to override.")
        sys.exit(1)
    # cannot use .apply on an empty dataframe without losing columns (which you need for "basin" filter below).
    if not df.empty:
        df = df.groupby(unique_track, as_index=False).apply(get_origmesh, initfile) # as_index=False avoids ValueError: 'basin' is both an index level and a column label, which is ambiguous. 
    # Keep if not a cyclogenesis track or TS strength for at least one time.
    df = df.groupby(unique_track).filter(lambda x: (x["basin"] != "TG").all() or (x["vmax"].max() >=  34))
    logging.info(f"after origmesh {df.groupby(unique_track).ngroups}")
    atcf.write(origmeshFile, df, debug=debug)
    fineprint += "\noutputfile: " + origmeshFile
    atcf_file = origmeshFile

fineprint += f"{warmcore_desc}\ninput file: {atcf_file}\ncreated {datetime.datetime.now()}" # outside plotbasin loop

# Regional plots
extents = atcf.basin_bounds
for plotbasin in plotbasins + ['storm']:
    plt.subplots(figsize=(11,7.5))
    plt.tight_layout(rect=(0,0.1,0.98,0.95)) # rect=(left,bottom,right,top) keeps room for legend and prevents labels from being chopped off
    projection = cartopy.crs.PlateCarree()
    # New axes for each plotbasin, or lat/lon grid spacing will be messed up.
    ax = atcf.get_ax(projection=projection)

    label_interval_hours = 24

    for (basin, cy, initial_time, model), track in df.groupby(unique_track):
        start_label = cy + "\n" + track.valid_time.min().strftime('%-m/%-d %-Hz') # used to have initial time. why?
        end_label = cy + " " + track.valid_time.max().strftime('%-m/%-d %-Hz')
        logging.info(f'plot_track {start_label}-{end_label}')
        atcf.plot_track(ax, start_label, track, end_label, label_interval_hours=label_interval_hours, debug=debug)
        # get storm name if present
        stormname = atcf.get_stormname(track)
        if stormname:
            # get ibtracks best track 
            idf, ifile = ibtracs.get_atcf(stormname, track.valid_time.min().strftime('%Y'), debug=debug)
            atcf.plot_track(ax, stormname, idf, stormname, label_interval_hours=label_interval_hours)

    if plotbasin in extents:
        extent = extents[plotbasin]
        ax.set_extent(extents[plotbasin])
        logging.info(f"plot basin {plotbasin} {extent}")

    # Set title
    ax.set_title(title)
    legax = atcf.TClegend(ax)

    plt.annotate(text=fineprint, xy=(10,2), xycoords='figure pixels', fontsize=5)

    ofile = atcf_file + f".{plotbasin}.png"
    plt.savefig(ofile)
    logging.info(f'created {os.path.realpath(ofile)}')
    #put on web server
    if args.to_server:
        cmd = [f"rsync {pngfile} ahijevyc@whitedwarf.mmm.ucar.edu:/web/htdocs/projects/mpas/plots/{init_yyyymmddhh}/{project}.gfdl_tracks.{basin}.png"]  
        logging.debug(cmd)
        subprocess.call(cmd)
    plt.clf()
