import argparse
import logging
import os
import pdb
import sys

import pandas as pd
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
parser.add_argument("ifile", help="input file")
parser.add_argument("-o", "--ofile", help="output file")
parser.add_argument("-d", "--debug", action="store_true")


args = parser.parse_args()
clobber = args.clobber

if args.ofile:
    ofile = args.ofile
else:
    ofile = args.ifile + ".new"
debug = args.debug

logging.debug(args)

assert clobber or not os.path.exists(
    ofile
), f"output file {ofile} exists. use --clobber option to override."


ds = xarray.open_dataset(args.ifile)

# unfortunately swap_dims loses the units so we can't use them when selecting
# below for tmean.
logging.warning("swap_dims")
swap_dims = dict(
    nIsoLevelsT="t_iso_levels",
    nIsoLevelsZ="z_iso_levels",
    nIsoLevelsU="u_iso_levels",
)
for k in swap_dims:
    if k in ds:
        ds = ds.swap_dims(k=swap_dims[k])

if "z_iso_levels" in ds:
    for level in ds.z_iso_levels:
        ilevel = str(int(level / 100))
        new_name = f"height_{ilevel}hPa"
        ds[new_name] = ds.z_isobaric.sel(z_iso_levels=level)

if "t_iso_levels" in ds:
    for level in ds.t_iso_levels:
        ilevel = str(int(level / 100))
        new_name = f"temperature_{ilevel}hPa"
        ds[new_name] = ds.t_isobaric.sel(t_iso_levels=level)

if "u_iso_levels" in ds:
    logging.warning("u_iso_levels")
    for level in ds.u_iso_levels:
        ilevel = str(int(level / 100))
        new_name = f"uzonal_{ilevel}hPa"
        ds[new_name] = ds.uzonal_isobaric.sel(u_iso_levels=level)
        new_name = f"umeridional_{ilevel}hPa"
        ds[new_name] = ds.umeridional_isobaric.sel(u_iso_levels=level)
    ds = ds.drop_vars("u_iso_levels")
else:
    logging.warning("u_iso_levels not in Dataset")


temperature_levels = [v for v in ds.variables if v.startswith("temperature") and v.endswith("hPa")]
tmean_levels = [v for v in temperature_levels if int(v[12:-3]) <= 500 and int(v[12:-3]) >= 300]
if "t_isobaric" in ds:
    # Instead of using meanT_500_300 from MPAS diagnostics just average the 5 levels with ncwa.
    # That is actually what hwrf_tave.exe does. (none of the dp-weighting used in MPAS diagnostics).
    logging.warning("tmean")
    ds["tmean"] = ds["t_isobaric"].sel(t_iso_levels=slice(30000,50000)).mean(dim="t_iso_levels")
elif tmean_levels:
    ds["tmean"] = ds[temperature_levels].to_dataarray().mean(dim="variable")
else:
    logging.warning("no tmean levels")
    

# change the way time is saved in netCDF file. Used to use cdo. cdo did some things like
# make time dimension UNLIMITED and add missing_value attribute to the variables, but I 
# don't think those things were necessary.
ds.time.encoding["units"] = f"hours since {ds.time.data[0]}"
# convert np.datetime64[ns] to pandas datetime. or else netcdf file has missing times.
ds["time"] = pd.to_datetime(ds.time)

# Tracker program gettrk_main.f v3.9a can't handle variables with a vertical dimension.
# It assumes a variable has only one pressure level.
# remove these variables.
drop_vars = ["z_isobaric","z_iso_levels","uzonal_isobaric","umeridional_isobaric", "t_isobaric"]
for var in drop_vars:
    if var in ds:
        ds = ds.drop_vars(var)

ds.to_netcdf(ofile)

print("wrote", ofile)
