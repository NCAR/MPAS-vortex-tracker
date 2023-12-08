"""
Read ATCF file
Optionally replace vitals with raw values from original mesh
Plot tracks meeting length and strength criteria.
Plots for different basins. 
"""
import argparse
import atcf
import cartopy.feature
import datetime
import ibtracs
import logging
import matplotlib.pyplot as plt
import mpas
import numpy as np
import pandas as pd
import pdb
import re, os
import subprocess
import sys

logger = logging.getLogger()
logger.setLevel(logging.INFO)

def main(args):

    debug = args.debug
    if debug:
        logger.setLevel(logging.DEBUG)

    logging.debug(args)

    alpha = args.alpha
    ax = args.ax
    basin       = args.basin
    clobber      = args.clobber
    diagdir      = args.diagdir
    extent       = args.extent
    initfile     = args.initfile
    initial_time = pd.to_datetime(args.initial_time)
    min_warmcore_percent = args.min_warmcore_percent
    models       = args.models
    origmesh     = args.origmesh
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
        assert not os.path.isfile(origmeshFile) or clobber, f"{origmeshFile} already exists. Use --clobber to override."

    df = atcf.read(atcf_file)

    # A track has a unique combination of basin, cy, initial_time, and model
    unique_track = ['basin', 'cy', 'initial_time', 'model']

    # Skip lines not from a specifically-requested initialization time.
    if initial_time:
        df = df[df.initial_time == initial_time]
        logging.info(f"kept {len(df)} {initial_time} initialization times")

    # Skip lines not from a specifically-requested model.
    if models:
        df = df[df['model'].isin(models)]
        logging.info(f"kept {len(df)} {models} lines.")


    # If you need original mesh values, initialize diagdir and initfile.
    if origmesh:
        if not diagdir:
            # Derive diagdir from atcf_file.
            m = re.search("/[12]\d{9}/", atcf_file)
            assert m, f"No diagdir supplied, and failed to derive from atcf_file: {atcf_file}"
            # path to diagnostic file original mesh
            diagdir = atcf_file[:m.end()]

        if not initfile:
            # Derive initfile from diagdir.
            # path to init.nc file with latCell,lonCell,etc.
            initfile = os.path.join(diagdir,"init.nc")

    warmcore_desc = f"warm core >= {min_warmcore_percent}% of the time"

    # If one of unique_tracks columns is constant, put in title. 
    title = ''
    for col in unique_track:
        if df[col].value_counts().count() == 1:
            # Add this column value to title. It doesn't vary.
            title += " " + col + "= " + str(df[col].iloc[0])

    # If warmcore column exists make sure at least one time is warmcore or unknown.
    if "warmcore" in df.columns:
        df = df.groupby(unique_track).filter(atcf.iswarmcore, min_warmcore_percent=min_warmcore_percent)
        logging.info(f'{df.groupby(unique_track).ngroups} after warmcore filter')
    else:
        warmcore_desc = f"No warm core filter"

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
            og, initfile = mpas.origmesh(track, initfile, diagdir, wind_radii_method=wind_radii_method)
        elif (model == 'ECMF' or re.compile("EE[0-9][0-9]$").match(model)):
            assert (track[["rad1","rad2","rad3","rad4"]] == 0).all().all(), (
                    "That is odd. You requested original mesh values but the atcf file already has non-zero radii")
            logging.info(f"{model} getting radius of max wind and radii of 34/50/64 kt wind from {diagdir}")
            og = atcf.origgrid(track, diagdir, ensemble_prefix=["PF","CF","EE","EE0"], wind_radii_method=wind_radii_method)
        elif "kowaleski/" in diagdir:
            # Tack grid nest number to end of original mesh output filename. 
            # Might help prevent applying d03 grid to d02 tracks, or vice versa.
            grid="d03"
            origmeshFile += grid
            og = atcf.origgridWRF(track, diagdir, grid=grid, wind_radii_method=wind_radii_method)
        else:
            logging.warning(f"no origmesh for {model}")
            return track # used to exit return code=3 but why do that if you have OFCL forecast in same file as model forecasts? that's fine.
        return og # replace track with og DataFrame so when track is plotted below, it is the original mesh stuff.

    if origmesh:
        # cannot use .apply on an empty dataframe without losing columns (which you need for "basin" filter below).
        if not df.empty:
            df = df.groupby(unique_track, as_index=False).apply(get_origmesh, initfile) # as_index=False avoids ValueError: 'basin' is both an index level and a column label, which is ambiguous. 
        # Keep if not a cyclogenesis track or TS strength for at least one time.
        df = df.groupby(unique_track).filter(lambda x: (x["basin"] != "TG").all() or (x["vmax"].max() >=  34))
        logging.info(f"after origmesh {df.groupby(unique_track).ngroups}")
        atcf.write(origmeshFile, df)
        atcf_file = origmeshFile


    # Regional plots

    if ax is None:
        fig = plt.figure(figsize=(11,7.5))
        ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax = atcf.decorate_ax(ax)

    if extent:
        ax.set_extent(extent)
    else:
        extent = atcf.basin_bounds[basin]
        ax.set_extent(extent)
    
    logging.info(f"plot extent {extent}")

    for (_, cy, initial_time, model), track in df.groupby(unique_track):
        if model == "OFCL":
            logging.warning("skip OFCL")
            continue
        start_label = None
        end_label = None
        duration  = track.valid_time.max() - track.valid_time.min()
        logging.info(f'plot_track {cy} {model} {duration}')
        atcf.plot_track(ax, start_label, track, end_label,
                label_interval_hours=None,
                alpha=alpha
                )
    # plot BEST track
    stormname = atcf.get_stormname(track)
    if stormname:
        logging.warning(f"get {stormname} best track") 
        idf, _ = ibtracs.get_df(stormname=stormname, year = track.valid_time.min().year)
        atcf.plot_track(ax, stormname, idf, stormname,
                label_interval_hours=24,
                scale=3)


    # Set title
    ax.set_title(title)

    if args.ax is None:
        finish(args)

    return ax

def get_parser():
    parser = argparse.ArgumentParser(description='plot tropical cyclone track/intensity from ATCF file',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--alpha', type=float, default=1, help='alpha channel (0-1)')
    parser.add_argument('atcf_file', type=argparse.FileType('r'), help='path to ATCF file to plot and possibily modify')
    parser.add_argument('--ax', help='plot in this axes object')
    parser.add_argument('--clobber', action="store_true", help='overwrite old file')
    parser.add_argument('-d','--debug', action="store_true", help='print debug messages')
    parser.add_argument('--diagdir', help='path to directory with original model diagnostic files. In case of ensemble, path to directory above member subdirectories i.e. do not include member subdirectory.')
    parser.add_argument('--initfile', help='path to init.nc file with latCell,lonCell,nEdgesOnCell,cellsOnCell')
    parser.add_argument('--initial_time', help='initialization date/hour')
    parser.add_argument('--min_warmcore_percent', type=float, default=25, help='minimum %% of warm core track times') # double percent sign to print single percent sign
    parser.add_argument('--models', nargs='*', help='only plot this model(s)')

    origmesh_group = parser.add_mutually_exclusive_group()
    origmesh_group.add_argument('--origmesh',    dest='origmesh', action="store_true", help='use original mesh')
    origmesh_group.add_argument('--no-origmesh', dest='origmesh', action="store_false", help='do not use original mesh')

    extent_group = parser.add_mutually_exclusive_group()
    extent_group.add_argument('--basin', help='basin to plot')
    extent_group.add_argument('-e', '--extent', nargs=4, type=float, help='minlon maxlon minlat maxlat (overrides basin)')

    parser.add_argument('--wind_radii_method', type=atcf.iswind_radii_method, default='max', help='method to get wind extent from original mesh')
    return parser

def finish(args):
    legax = atcf.TClegend()
    fineprint = f"{args} created {datetime.datetime.now()}" # outside basin loop
    plt.annotate(text=fineprint, xy=(10,2), xycoords='figure pixels', fontsize=5, wrap=True)
    atcf_file = os.path.realpath(args.atcf_file.name)
    sfx = f".{args.basin}.png" if args.basin else ".png"
    ofile = atcf_file + sfx
    plt.savefig(ofile, dpi=120)
    logging.warning(ofile)


if __name__ == "__main__":
    # tried to put in main() but run_plot_atcf.py calls main(args) to pass ifile, alpha, extent, ax, etc...
    parser = get_parser()
    args = parser.parse_args()
    main(args)
