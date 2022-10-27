import argparse
import logging
import os
import pdb
import sys
import xarray

# =============Arguments===================
parser = argparse.ArgumentParser(description = "Unstack variables with a vertical dimension into multiple single-level variables with vertical level encoded in the variable name", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--clobber", action='store_true', help='clobber output file even if it exists')
parser.add_argument("ifile", type=argparse.FileType('r'), help='input file')
parser.add_argument("-o", "--ofile", type=argparse.FileType('w'), help='output file')
parser.add_argument("-d", "--debug", action='store_true')


args = parser.parse_args()
clobber     = args.clobber
ifile       = os.path.realpath(args.ifile.name)
if args.ofile:
    ofile       = os.path.realpath(args.ofile.name)
else:
    ofile = ifile + ".new"
debug       = args.debug

if debug:
    ofile = ofile + "_debug"
    print(args)

if os.path.exists(ofile) and not clobber:
    print("output file",ofile, "exists. Exiting.")
    print("use --clobber option to override.")
    sys.exit(1)


ds = xarray.open_dataset(ifile)
ds = ds.swap_dims(dict(nIsoLevelsT="t_iso_levels",nIsoLevelsZ="z_iso_levels"))

for level in ds.z_iso_levels:
    ilevel = str(int(level/100))
    new_name = f"height_{ilevel}hPa"
    ds[new_name] = ds.z_isobaric.sel(z_iso_levels=level)
ds = ds.drop_vars("z_iso_levels")

for level in ds.t_iso_levels:
    ilevel = str(int(level/100))
    new_name = f"temperature_{ilevel}hPa"
    ds[new_name] = ds.t_isobaric.sel(t_iso_levels=level)
#ds = ds.drop_vars("t_iso_levels") # needed for average temperature derivation later

if "u_iso_levels" in ds:
    for level in ds.u_iso_levels:
        ilevel = str(int(level/100))
        new_name = f"uzonal_{ilevel}hPa"
        ds[new_name] = ds.uzonal_isobaric.sel(u_iso_levels=level)
        new_name = f"umeridional_{ilevel}hPa"
        ds[new_name] = ds.umeridional_isobaric.sel(u_iso_levels=level)
    ds = ds.drop_vars("u_iso_levels")
else:
    logging.warning("u_iso_levels not in Dataset")


# uncomment if you want to remove these variables.
#ds = ds.drop(["z_isobaric","uzonal_isobaric","umeridional_isobaric", "t_isobaric"])
ds.to_netcdf(ofile)
print("wrote",ofile)
