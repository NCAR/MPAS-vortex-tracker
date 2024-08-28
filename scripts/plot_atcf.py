"""
Read ATCF file
Optionally replace vitals with raw values from original mesh
Plots for different basins. 
"""
import argparse
import cartopy
import os
import sys
from atcf import (
        basin_bounds,
        decorate_ax,
        iswind_radii_method,
        origgrid,
        plot_track,
        read,
        write,
        TClegend,
        )
import logging
import matplotlib.pyplot as plt
import pandas as pd
import pdb
import re

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", force=True)

def bestTrack(stormname:str) -> pd.DataFrame:
    from stormevents.nhc import VortexTrack
    # BEST track
    tfile = os.path.join(os.getenv("TMPDIR"), f"{stormname}.csv")
    if os.path.exists(tfile):
        logging.warning(f"read {tfile}")
        idf = pd.read_csv(tfile, parse_dates=["valid_time","track_start_time"], index_col=0)
    else:
        logging.warning(f"get {stormname} best track") 
        idf = VortexTrack(storm=stormname, advisories=["BEST"]).data.rename(
                columns={
                    "datetime" : "valid_time",
                    "latitude" : "lat",
                    "longitude" : "lon",
                    "max_sustained_wind_speed" : "vmax",
                }
            )
        logging.warning(f"save {tfile}")
        idf.to_csv(tfile)
    # VortexTrack in knots
    return idf

def main(args):

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s", force=True)

    logging.debug(args)

    ax = args.ax
    diagdir = args.diagdir
    extent = args.extent
    initfile = args.initfile
    origmesh = args.origmesh
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
        assert not os.path.isfile(origmeshFile) or args.clobber, f"{origmeshFile} already exists. Use --clobber to override."

    df = read(atcf_file)

    # A track has a unique combination of basin, cy, initial_time, and model
    unique_track = ['basin', 'cy', 'initial_time', 'model']

    # If you need original mesh values, initialize diagdir and initfile.
    if origmesh:
        if not diagdir:
            # Derive diagdir from atcf_file.
            m = re.search(r"/[12]\d{9}/", atcf_file)
            assert m, f"No diagdir supplied, and failed to derive from atcf_file: {atcf_file}"
            # path to diagnostic file original mesh
            diagdir = atcf_file[:m.end()]

        if not initfile:
            # Derive initfile from diagdir.
            # path to init.nc file with latCell,lonCell,etc.
            initfile = os.path.join(diagdir,"init.nc")


    # If one of unique_tracks columns is constant, put in title. 
    title = ''
    for col in unique_track:
        if df[col].value_counts().count() == 1:
            # Add this column value to title. It doesn't vary.
            title += " " + col + "= " + str(df[col].iloc[0])

    def get_origmesh(track, initfile):
        basin, cy, initial_time, model = track.name
        # Do you need to find the original mesh values? Not if track has defined originalmeshfile . 
        if 'originalmeshfile' in track.columns and track["originalmeshfile"].all():
            logging.warning(f"{track_id} already has original mesh values.")
            return track
        elif (model == 'ECMF' or re.compile("EE[0-9][0-9]$").match(model)):
            assert (track[["rad1","rad2","rad3","rad4"]] == 0).all().all(), (
                    "That is odd. You requested original mesh values but the atcf file already has non-zero radii")
            logging.info(f"{model} getting radius of max wind and radii of 34/50/64 kt wind from {diagdir}")
            og = origgrid(track, diagdir, ensemble_prefix=["PF","CF","EE","EE0"], wind_radii_method=wind_radii_method)
        elif (model == 'MPAS'):
            import mpas
            og, initfile = mpas.origmesh(track, initfile, diagdir, wind_radii_method=wind_radii_method)
        else:
            logging.warning(f"no origmesh for {model}")
            return track # used to exit return code=3 but why do that if you have OFCL forecast in same file as model forecasts? that's fine.
        return og # replace track with og DataFrame so when track is plotted below, it is the original mesh stuff.

    if origmesh:
        # cannot use .apply on an empty dataframe without losing columns (which you need for "basin" filter below).
        if not df.empty:
            df = df.groupby(unique_track).apply(get_origmesh, initfile, include_groups=False)
        df = df.reset_index(drop=False)
        # Keep if not a cyclogenesis track or TS strength for at least one time.
        df = df.groupby(unique_track).filter(lambda x: (x.name[unique_track.index("basin")] != "TG") or (x["vmax"].max() >=  34))
        logging.info(f"after origmesh {df.groupby(unique_track).ngroups}")
        write(df, origmeshFile)
        atcf_file = origmeshFile


    # Regional plots

    if ax is None:
        logging.warning("ax is None, so create fig and ax")
        fig = plt.figure(figsize=(10,7.5))
        ax = plt.axes(projection=cartopy.crs.PlateCarree())
    ax = decorate_ax(ax)

    if extent:
        ax.set_extent(extent)
    else:
        extent = basin_bounds[args.basin]
        ax.set_extent(extent)
    
    logging.info(f"plot extent {extent}")

    # Convert ECMWF ensemble Vmax from 10-min to 1-min average
    ten2one = atcf_file.endswith("scaled_vmax") or atcf_file.endswith("scaled_vmax.Knaff_Zehr_pmin")
    for (basin, cy, initial_time, model), track in df.groupby(unique_track):
        linewidth = 2
        linestyle = "solid"
        alpha = args.alpha
        if model == "OFCL":
            linewidth = 3
            linestyle = "dashed"
            alpha = 1
            ten2one = False # Assume OFCL Vmaxl is already 1-min average
        if ten2one:
            import nos
            logging.warning(f"convert {model} vmax 10-min average to 1-min") 
            track["vmax"] = track["vmax"] / nos.one2ten
        start_label = None
        end_label = None
        duration  = track.valid_time.max() - track.valid_time.min()
        logging.info(f'plot_track {cy} {model} {duration}')
        plot_track(ax, start_label, track, end_label,
                label_interval_hours=24,
                alpha=alpha,
                linewidth=linewidth,
                linestyle=linestyle,
                )
    if args.stormname:
        idf = bestTrack(args.stormname)
        plot_track(
                    ax,
                    "",
                    idf, 
                    "",
                    label_interval_hours=24,
                    linewidth = 3,
                    linestyle = "solid",
                    scale=2,
                    alpha=0.4,
                )


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
    parser.add_argument('--ofile', help='output filename')
    parser.add_argument('--stormname', help='plot best track. name like al02024')

    origmesh_group = parser.add_mutually_exclusive_group()
    origmesh_group.add_argument('--origmesh',    dest='origmesh', action="store_true", help='use original mesh')
    origmesh_group.add_argument('--no-origmesh', dest='origmesh', action="store_false", help='do not use original mesh')

    extent_group = parser.add_mutually_exclusive_group()
    extent_group.add_argument('--basin', help='basin to plot')
    extent_group.add_argument('-e', '--extent', nargs=4, type=float, help='minlon maxlon minlat maxlat (overrides basin)')

    parser.add_argument('--wind_radii_method', type=iswind_radii_method, default='max', help='method to get wind extent from original mesh')
    return parser

def finish(args):
    """add legend and fine print"""
    legax = TClegend(aspect=30, pad=0.08)
    fineprint = f"{args} created {pd.Timestamp.now()}" 
    plt.annotate(text=fineprint, xy=(10,2), xycoords='figure pixels', fontsize=5, wrap=True)
    atcf_file = os.path.realpath(args.atcf_file.name)
    sfx = f".{args.basin}.png" if args.basin else ".png"
    if args.ofile:
        ofile = args.ofile
    else:
        ofile = atcf_file + sfx
    plt.savefig(ofile, dpi=120)
    logging.warning(ofile)


if __name__ == "__main__":
    # tried to put in main() but run_plot_atcf.py calls main(args) to pass ifile, alpha, extent, ax, etc...
    logging.warning("get parser")
    parser = get_parser()
    args = parser.parse_args()
    logging.warning("start main()")
    main(args)
