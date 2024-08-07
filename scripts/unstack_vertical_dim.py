import argparse
import logging
import os
import pdb
import sys

import metpy
from metpy.units import units
import xarray

# =============Arguments===================
parser = argparse.ArgumentParser(
    description="Unstack variables with a vertical dimension into multiple single-level "
    "variables with vertical level encoded in the variable name",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "--clobber", action="store_true", help="clobber output file even if it exists"
)
parser.add_argument("ifile", type=argparse.FileType("r"), help="input file")
parser.add_argument("-o", "--ofile", type=argparse.FileType("w"), help="output file")
parser.add_argument("-d", "--debug", action="store_true")


args = parser.parse_args()
clobber = args.clobber
ifile = os.path.realpath(args.ifile.name)
if args.ofile:
    ofile = os.path.realpath(args.ofile.name)
else:
    ofile = ifile + ".new"
debug = args.debug

if debug:
    ofile = ofile + "_debug"
    print(args)

assert clobber or not os.path.exists(
    ofile
), f"output file {ofile} exists. use --clobber option to override."


ds = xarray.open_dataset(ifile).metpy.quantify()

# unfortunately swap_dims loses the units so we can't use them when selecting
# below for tmean.
ds = ds.swap_dims(
    dict(
        nIsoLevelsT="t_iso_levels",
        nIsoLevelsZ="z_iso_levels",
        nIsoLevelsU="u_iso_levels",
    )
)

for level in ds.z_iso_levels:
    ilevel = str(int(level / 100))
    new_name = f"height_{ilevel}hPa"
    ds[new_name] = ds.z_isobaric.sel(z_iso_levels=level)

for level in ds.t_iso_levels:
    ilevel = str(int(level / 100))
    new_name = f"temperature_{ilevel}hPa"
    ds[new_name] = ds.t_isobaric.sel(t_iso_levels=level)

if "u_iso_levels" in ds:
    for level in ds.u_iso_levels:
        ilevel = str(int(level / 100))
        new_name = f"uzonal_{ilevel}hPa"
        ds[new_name] = ds.uzonal_isobaric.sel(u_iso_levels=level)
        new_name = f"umeridional_{ilevel}hPa"
        ds[new_name] = ds.umeridional_isobaric.sel(u_iso_levels=level)
    ds = ds.drop_vars("u_iso_levels")
else:
    logging.warning("u_iso_levels not in Dataset")


# Instead of using meanT_500_300 from MPAS diagnostics just average the 5 levels with ncwa.
# That is actually what hwrf_tave.exe does. (none of the dp-weighting used in MPAS diagnostics).
ds["tmean"] = ds["t_isobaric"].sel(t_iso_levels=slice(30000,50000)).mean(dim="t_iso_levels")

# change the way time is saved in netCDF file. Used to use cdo. cdo did some things like
# make time dimension UNLIMITED and add missing_value attribute to the variables, but I 
# don't think those things were necessary.
ds.time.encoding["units"] = f"hours since {ds.time.data[0]}"

# Tracker program gettrk_main.f v3.9a can't handle variables with a vertical dimension.
# It assumes a variable has only one pressure level.
# remove these variables.
ds = ds.drop_vars(["z_isobaric","z_iso_levels","uzonal_isobaric","umeridional_isobaric", "t_isobaric"])

ds.metpy.dequantify().to_netcdf(ofile)
print("wrote", ofile)
