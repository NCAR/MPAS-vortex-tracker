# MPAS-vortex-tracker
track tropical cyclones in MPAS

## Preprocessing

### Smooth and interpolate to lat-lon grid

Clone github repository: `git clone https://github.com/NCAR/MPAS-vortex-tracker.git`

datetime-fortran module came from [https://github.com/wavebitscientific/datetime-fortran]. 

To compile:
```
module reset
module load gcc
cd MPAS-vortex-tracker/src
make clean
make
cd ..
```

creates executable: `mpas_to_latlon` 

If you have never interpolated a particular mesh before, `mpas_to_latlon` will first create a save file.  The save file will contain expensive triangulation information and a time variable. To create the save file, run `mpas_to_latlon` on a file that has the variables latCell, lonCell and xtime. init.nc is such a file. `mpas_to_latlon` will create the save file in the path defined in the environmental variable TMPDIR. Set the environmental variable `$TMPDIR` and make it with mkdir, if it doesn't exist already. Finally, execute `mpas_to_latlon`. Ask for a field that is in `init.nc`, like precipw. 

```
setenv TMPDIR /glade/scratch/$USER/temp
mkdir -p $TMPDIR
echo precipw | bin/mpas_to_latlon init.nc $TMPDIR/latlon.init.nc 0.5 25 mpas3 -5 50 0
```

The echo command pipes the field name `precipw` into `mpas_to_latlon`. In this case, `mpas_to_latlon` will create two files. The first file is a save file in `$TMPDIR`. It is called `mpas3_02621442_0.500deg -5.000N+50.000N  +0.000E`. The min and max latitude and the min longitude portions of the file name are fixed-width format, so they may contain spaces. For example, in this case, there is a space before the southernmost latitude portion `-5.000N` of the file name and two spaces before the min longitude `lon0` portion of the file name. `02621442` is the number of cells in the mesh (zero-padded). The save file will be reused next time this `mpas3` mesh is interpolated. The second file is the interpolated output, `$TMPDIR/latlon.init.nc`. It may be deleted, but could be useful for debugging. 

To interpolate multiple fields, provide a list of fields to mpas_to_latlon.  These may be typed in manually or be in a text file, one to a line, like this:

```
t_isobaric
meanT_500_300
mslp
u10
umeridional_isobaric
uzonal_isobaric
v10
z_isobaric
```

#### Usage
`cat scripts/mpas_fields_to_interpolate.txt | grep -v # | bin/mpas_to_latlon fname outfile grid_spacing filter_radius_km meshid [lat0 [lat1 [startLon]]]`

```
 fname            : input MPAS filename
 outfile          : output filename
 grid_spacing     : output grid spacing in degrees
 filter_radius_km : radius of circular smoothing filter in km
 meshid           : string to identify MPAS mesh (uni, wp, us, etc.)
 lat0             : minimum latitude in deg (default -90)
 lat1             : maximum latitude in deg (default +90)
 startLon         : starting longitude for output grid (default -180)
```

#### Example:
```cat mpas_fields_to_interpolate.txt | mpas_to_latlon diagnostics.2013-09-01_00:00:00.nc latlon_0.500deg_25km/outfile.nc 0.5 25 mpas3 -5 50 0```

This example interpolates the list of fields in mpas_fields_to_interpolate.txt to a 0.5 degree lat-lon grid after smoothing with a radius of influence of 25 km.  meshid is the arbitrary string label describing the MPAS mesh. Choosing a unique string can help distinguish different basins, even if the meshes have the same number of cells and the same lat-lon output parameters.  For example, I use 'wp' for the western Pacific mesh and 'ep' for the eastern Pacific mesh. In this example, the string description of the MPAS mesh is mpas3 and the latitude range is -5 to +50 N. The output file is called output.nc, but in practice, it is helpful to use a more descriptive name or subdirectory like latlon_0.500deg_25km/diag.2013-09-01_00:00:00.nc. Keep in mind, if the subdirectory does not exist already, or you don't have permission to write to the subdirectory, mpas_to_latlon will exit with an error.

Smoothing is performed iteratively on the native MPAS mesh.  During each pass the cells are averaged with all cell neighbors that share an edge.  Each cell is assigned a smoothing pass count (fpasses) that is inversely proportional to its radius. Small cells get averaged with their neighbors more times than large cells. This is like smoothing with a spatial filter of constant size. 

Smoothing algorithm pseudo code:

```
! Input : field(nCells) (unsmoothed values on cells)
!       : areaCell(nCells) (area of cell in sq. km)
!       : roi (radius of influence in km)
! Output: field(nCells) (smoothed values)
! 
! fpasses(iCell): number of smoothing passes for each cell
do iCell = 1, nCells
   fpasses(iCell) = roi/sqrt(areaCell(iCell)/pi)
end do

do fpass = 1, max(fpasses)
   do iCell = 1, nCells
      tmp(iCell) = field(iCell)
      ! new cell value is half old value and
      ! half average of neighbors.
      tmp(iCell) = 0.5 * field(icell) + 0.5 * avg(neighbors)
   end do
   do iCell = 1, nCells
	field(iCell) = tmp(iCell)
   end do
end do
```

It's a linear interpolation from the mpas cell centers to the lat-lon grid cell centers. Delaunay triangulation is used to find the closest 3 mesh cell centers and obtain weights for the final interpolated values.  For vorticity, which is not defined on vertices, the interpolation is preceded by a step where the values on the cell vertices are averaged to get cell-centered values. Then the smoothing and interpolation proceed as normal. 

A shell script scripts/run_mpas_to_latlon.csh runs mpas_to_latlon on all the diagnostic files for one model run. 

#### Usage

`scripts/run_mpas_to_latlon.csh [-i idir] [-w working_dir] [--mesh mesh_id] [--delta d] [--executable EXECUTABLE] [-t fields_to_interpolate] [--lat0 lat0] [--lat1 lat1] [--filter_radius_km filter_radius_km] [--lon0 lon0]`

#### Defaults:
```
working_dir                      : $TMPDIR
mesh_id                          : tk707_conus
d (grid spacing in degrees)      : 0.5 (0.5°)
EXECUTABLE                       : $SCRATCH/MPAS-vortex-tracker/bin/mpas_to_latlon
fields_to_interpolate file       : $SCRATCH/MPAS-vortex-tracker/scripts/mpas_fields_to_interpolate.txt
southern latitude (lat0)         : -5 (5°S)
northern latitude (lat1)         : 55 (55°N)
starting longitude (lon0)        : 0 (0°E)
filter radius (filter_radius_km) : 25 (25 km)
Arguments may be in any order.
```

#### Example: 
```scripts/run_mpas_to_latlon.csh --mesh mpas15```

This will Interpolate all `diag*.nc` files in `$TMPDIR/mpas15/` to a lat-lon grid. 

## Convert to format that vortex tracker can read

As of version 3.9.1, the vortex tracker can process netCDF input.  Longitude must be 0-360, not -180 to +180 or else the longitude is output in atcf with the 'E' suffix, even with negative values. This is handled by the script `scripts/mpas_ll_gettrk.csh`.  The functions of `mpas_ll_gettrk.csh` are described below.

`mpas_ll_gettrk.csh`

  - Change to the directory with interpolated MPAS lat/lon diagnostic files.
  
  - Make GFDL tracker working directories. ($trackertype can be tracker or tcgen)
  `mkdir -p gfdl_tracker/$trackertype`
  
  - Create fort.15 file. This is a numbered list of forecast lead times in minutes (index of file starting with "0001" and forecast lead time in minutes). Use forecast_hour_links.py. Handles missing diag files. 
  
  - Concatenate the MPAS lat/lon diagnostic files into a single file:
  `ncrcat ../diag*.nc all.nc`
  
  - Unstack wind, height, and temperature variables. GFDL tracker doesn't expect variables with a vertical dimension.
  `python scripts/unstack_vertical_dim.py all.nc --ofile diag.$meshid.$ymd$h.nc --clobber`
  
  - If you run the tracker in ‘tcgen’ (tropical cyclone genesis) mode, you can skip the rest of the section. Otherwise for ‘tracker’ mode, you need the initial storm positions in a file called tcvit_rsmc_storms.txt.  wget_tcvitals.csh will download TC vitals for the requested init time and put it in the aforementioned file.
  `scripts/wget_tcvitals.csh $ymd $h tcvit_rsmc_storms.txt`
  where `$ymd` is the initialization time in yyyymmdd format and `$h` is hour. 
  
  - Make empty fort.14 file.
  `touch fort.14`
  
  - Fill namelist.
  
  - Run tracker.
  `$BINDIR/gettrk.exe < namelist > log`

# Run tracker on GFS - <a href="https://docs.google.com/document/d/1vgNUB4GW0FOgpD3tZUm8GQ2jV8ECW5SI7tfZVcX12UU/edit?usp=sharing">Google Doc</a>

## Matching

## Stats
