#!/bin/csh
#
# convert MPAS lat/lon files to input for GFDL vortex tracker
# remember to interpolate to lat/lon with run_mpas_to_latlon.csh first

# usage
# mpas_ll_gettrk.csh [-w model] [-t (tracker|tcgen)] [-d] [-i workdir_parent] [--justplot]

module load conda
conda activate npl-2024b
module load nco
setenv BINDIR $SCRATCH/standalone_gfdl-vortextracker/trk_exec
setenv EXEDIR $SCRATCH/MPAS-vortex-tracker


# set defaults
set debug=0
set dxdetails=_0.500deg_025km
set justplot=0
set meshid=tk707_conus
set trackertype=tcgen # trackertype=tracker is traditional and only tracks stuff mentioned at t=0.  tcgen does genesis too
set workdir=$TMPDIR

while ("$1" != "")
	if ("$1" == "--dx") then # optional argument --dx grid spacing and smoothing directory. 
		shift
		set dxdetails="$1"
	endif
	if ("$1" == "--mesh") then 
		shift
		set meshid="$1"
	endif
	if ("$1" == "-t") then # optional argument -t can determine trackertype. (tracker or tcgen)
		shift
		set trackertype="$1"
	endif
	if ("$1" == "-d") set debug=1
	if ("$1" == "--justplot") set justplot=1
	if ("$1" == "-w") then # optional argument -w determines work directory. 
		shift
		set workdir="$1"
	endif
	shift
end

if ($debug) set echo

if ("$1" != "") then
	echo unknown argument $1
	exit
endif

# make north and south boundaries just inside boundaries of diag.latlon files.
source ./bounds.csh
set northbd=`echo "$latmax - 5"|bc`
set southbd=`echo "$latmin + 5"|bc`

cd $workdir/$meshid/latlon$dxdetails
if ($justplot) goto PLOT


mkdir -p gfdl_tracker/$trackertype
# used to use simple wildcard, but it matched hours between multiples of 3
# Hourly files were present for random dates, like mpas3/2013092700.
set fort15=gfdl_tracker/fort.15
ls diag*.nc|python $EXEDIR/scripts/forecast_hour_links.py > $fort15

# Enter directory for vortex tracker files
cd gfdl_tracker



# Concatenate forecast lead times into 1 file. Don't reuse the name all.nc or treat as temporary file.
# Read later for mslp.
ncrcat -O ../diag*.nc all.nc

# Used to assume workdir contained a yyyymmddhh string and pull initdate from that. Now get initdate from earliest 
# time in all.nc using ncdump. 
# Find last line with " time = " in it, grab its first field split by commas, grab its 2nd field split 
# by equals sign, and remove quotes. Extract time components with date function and various output formats.
set d=`ncdump -v time all.nc -t | grep " time = " | tail -n 1 | cut -d , -f 1 | cut -d= -f2 | tr -d \"`
set initdate=`date --date="$d" "+%Y-%m-%d"`
set bcc=`echo $d|cut -c1-2`
set byy=`date --date="$d" "+%y"`
set m=`date --date="$d" "+%m"`
set dd=`date --date="$d" "+%d"`
set h=`date --date="$d" "+%H"`
set ymd=$bcc$byy$m$dd

set out=diag.$meshid.$ymd$h.nc

if (! -s $out) then
    echo Making $out
    # separate variables with a vertical dimension into multiple variables named by level.
    python $EXEDIR/scripts/unstack_vertical_dim.py all.nc --ofile $out --clobber
    if (! $debug) rm -v all.nc
endif



# Download combined TC vitals for requested year, save in 
# tcvit_rsmc_storms.txt (tracker, formerly fort.12).
# or tcvit_genesis_storms.txt (tcgen, formerly fort.12)- couldn't get this to work with tcvit_rsmc_storms.txt existing.
set tcvitalsfile=tcvit_rsmc_storms.txt
#if ($trackertype == 'tcgen') set tcvitalsfile=tcvit_genesis_storms.txt
# Extract the date and hour that match the current init time
# Consider using --fiorino option. 
$EXEDIR/scripts/wget_tcvitals.csh $ymd $h $tcvitalsfile
if($trackertype == 'tracker' && ! -s $tcvitalsfile )then
    echo "No storms at initial time $ymd$h. Exiting $0." > no_tracker_storms
    cp no_tracker_storms $trackertype
    continue
endif  

# Make empty fort.14 file.
touch fort.14

# Don't use 850 and 700 zeta. gettrk doesn't seem to derive them from u and v correctly, leaving zeta missing. 
set phs=vtt
cat <<NL > namelist

&datein
    inp%bcc=$bcc, inp%byy=$byy, inp%bmm=$m,
    inp%bdd=$dd, inp%bhh=$h, inp%model=1,
    inp%modtyp='global',
    inp%lt_units='hours',
    inp%file_seq='onebig',
    inp%nesttyp='fixed' 
/
&atcfinfo
      atcfnum=0, atcfname='MPAS',
      atcfymdh=$bcc$byy$m${dd}$h, atcffreq=600
/
&trackerinfo trkrinfo%westbd=0,
  trkrinfo%eastbd=360,
  trkrinfo%northbd=$northbd,
  trkrinfo%southbd=$southbd,
  trkrinfo%type='$trackertype',
  trkrinfo%mslpthresh=0.0015,
  trkrinfo%use_backup_mslp_grad_check='y',
  trkrinfo%v850thresh=1.5000,
  trkrinfo%use_backup_850_vt_check='y',
  trkrinfo%gridtype='global',
  trkrinfo%enable_timing=0,
  trkrinfo%contint=100.0,
  trkrinfo%want_oci=.FALSE.,
  trkrinfo%out_vit='y',
  trkrinfo%use_land_mask='n',
  trkrinfo%inp_data_type='netcdf'
/
&phaseinfo phaseflag='y',
       phasescheme='$phs',
       wcore_depth=1.0
/
&structinfo structflag='n',
        ikeflag='n'
/
&fnameinfo  gmodname='diag',
        rundescr='$meshid',
        atcfdescr=''
/
&waitinfo use_waitfor='n'
/
&netcdflist netcdfinfo%num_netcdf_vars=999,
  netcdfinfo%netcdf_filename='$out',
  netcdfinfo%rv850name='X', 
  netcdfinfo%rv700name='X',
  netcdfinfo%u850name='uzonal_850hPa',
  netcdfinfo%v850name='umeridional_850hPa',
  netcdfinfo%u700name='uzonal_700hPa',
  netcdfinfo%v700name='umeridional_700hPa',
  netcdfinfo%z850name='height_850hPa',
  netcdfinfo%z700name='height_700hPa',
  netcdfinfo%mslpname='mslp',
  netcdfinfo%usfcname='u10',
  netcdfinfo%vsfcname='v10',
  netcdfinfo%u500name='uzonal_500hPa',
  netcdfinfo%v500name='umeridional_500hPa',
  netcdfinfo%tmean_300_500_name='tmean',
  netcdfinfo%z500name='height_500hPa',
  netcdfinfo%z200name='height_200hPa',
  netcdfinfo%lmaskname='X',
  netcdfinfo%z900name='X',
  netcdfinfo%z800name='height_800hPa',
  netcdfinfo%z750name='height_750hPa',
  netcdfinfo%z650name='X',
  netcdfinfo%z600name='height_600hPa',
  netcdfinfo%z550name='height_550hPa',
  netcdfinfo%z450name='height_450hPa',
  netcdfinfo%z400name='height_400hPa',
  netcdfinfo%z350name='height_350hPa',
  netcdfinfo%z300name='height_300hPa',
  netcdfinfo%time_name='time',
  netcdfinfo%lon_name='lon',
  netcdfinfo%lat_name='lat',
  netcdfinfo%time_units='hours'
/
&parmpreflist user_wants_to_track_zeta850='n',
  user_wants_to_track_zeta700='n',
  user_wants_to_track_wcirc850='y',
  user_wants_to_track_wcirc700='y',
  user_wants_to_track_gph850='y',
  user_wants_to_track_gph700='y',
  user_wants_to_track_mslp='y',
  user_wants_to_track_wcircsfc='y',
  user_wants_to_track_zetasfc='y',
  user_wants_to_track_thick500850='y',
  user_wants_to_track_thick200500='y',
  user_wants_to_track_thick200850='y'
/
&verbose
verb=1
/
NL

echo "Running tracker..."
echo $BINDIR/gettrk.exe \< namelist \> log
$BINDIR/gettrk.exe < namelist > log
if ($status != 0) then
    tail log
    exit 1
endif

cp fort.15 fort.61 fort.62 fort.64 fort.68 fort.69 namelist log $tcvitalsfile $trackertype
if($trackertype == 'tcgen') mv fort.66 fort.67 fort.77 $trackertype

# jump here if --justplot option 
PLOT:
# Between tracker and tcgen, we prefer tcgen online.
set atcf_with_warmcore_column=fort.64

# Read tracks and plot with script
# -f force overwrite 
# The 1st basin sets the lat/lon grid tick interval for all basins.
set echo
# TODO: allow cp basin. when I added cp to the basin list, it collapsed the global plot longitude range to zero.
python ~ahijevyc/bin/plot_atcf.py $atcf_with_warmcore_column --basin global al wp ep io \
    --origmesh --diagdir ../ --initfile $meshid/init.nc --force_new 

