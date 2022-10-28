#!/bin/csh
#
# convert MPAS lat/lon files to input for GFDL vortex tracker
# remember to interpolate to lat/lon with run_mpas_to_latlon.csh first

# usage
# mpas_ll_gettrk.csh [yyyymmddhh] [-w model] [-t (tracker|tcgen)] [-d] [-i workdir_parent] [--justplot]

module load cdo
module load nco
module load python # for forecast_hour_links.py 
module load ncarenv # for ncar_pylib

printenv NCAR_HOST
if ( $NCAR_HOST =~ cheyenne* ) then
    set BINDIR = /glade/u/home/ahijevyc/bin_cheyenne
    ncar_pylib /glade/work/ahijevyc/20201220_cheyenne_daa # for cartopy and atcf modules used in plot_atcf.py as user other than ahijevyc, like mpasrt
endif 
if ( $NCAR_HOST =~ dav ) then
    set BINDIR = /glade/u/home/ahijevyc/bin_dav
    ncar_pylib /glade/work/ahijevyc/20201220_daa_casper # for cartopy and atcf modules used in plot_atcf.py as user other than ahijevyc, like mpasrt
endif 


setenv BINDIR /glade/scratch/$USER/standalone_gfdl-vortextracker_v3.9a/trk_exec
setenv EXEDIR /glade/scratch/$USER/MPAS-vortex-tracker


# set defaults
set dxdetails=_0.500deg_025km
set debug=0
set idir=/glade/scratch/mpasrt
set justplot=0
set trackertype=tcgen # tracker is traditional and only tracks stuff mentioned at t=0.  tcgen does genesis too
set yyyymmddhh=`date -u +%Y%m%d`

while ("$1" != "")
	if ("$1" =~ 20[0-9][0-9][01][0-9][0-3]*) set yyyymmddhh="$1"
	if ("$1" == "--dx") then # optional argument --dx grid spacing and smoothing directory. 
		shift
		set dxdetails="$1"
	endif
	if ("$1" == "-i") then # optional argument -i determines parent working directory. 
		shift
		set idir="$1"
	endif
	if ("$1" == "-w") then # argument -w determines working directory under $idir. added this because the # of potential working directory patterns is getting out of hand.
		shift
		set mp="$1"
	endif
	if ("$1" == "-t") then # optional argument -t can determine trackertype. (tracker or tcgen)
		shift
		set trackertype="$1"
	endif
	if ("$1" == "-d") set debug=1
	if ("$1" == "--justplot") set justplot=1
	shift
end

if ($debug) set echo

if ("$1" != "") then
	echo unknown argument $1
	exit
endif

umask 002

set workdir = $idir/$mp/$yyyymmddhh/latlon$dxdetails
if ($justplot) goto PLOT


cd $workdir

mkdir -p gfdl_tracker/$trackertype
# used to use simple wildcard, but it matched hours between multiples of 3
# Hourly files were present for random dates, like mpas3/2013092700.
set if=0
set fort15=gfdl_tracker/fort.15
if (-e $fort15) rm $fort15
foreach f (`ls diag.*$dxdetails.nc|sort` )

    # Append line to fort.15 (index of file starting with "0001" and forecast lead time in minutes) 
    # Used to be at end of foreach block but it needs to occur with every iteration of the diagnostics
    # file loop, before any "continue" clause.
    # Use forecast_hour_links.py to output the number of minutes. Allows missing diag files. 
    set fmin=`$EXEDIR/scripts/forecast_hour_links.py -m $f`
    printf '%04d %05d\n' `expr $if + 1` $fmin >> $fort15
    @ if++

end

# Enter directory for vortex tracker files
cd gfdl_tracker


# put lock file in working directory of tracker. 
# Don't want more than 1 tracker running in the same directory.
set lock=.mpas_ll_gettrk.lock
if (-e $lock) then
    echo found lock file $lock in `pwd` 
    exit 1
endif
touch $lock

# Concatenate forecast lead times into 1 file. Don't reuse the name all.nc or treat as temporary file.
# Read later for mslp.
if (! -s all.nc) ncrcat -O ../diag.*$dxdetails.nc all.nc

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

set out=diag.$mp.$ymd$h
if (! -s $out) then
    echo Making $out
    # separate levels
    python $EXEDIR/scripts/unstack_vertical_dim.py all.nc --ofile tmp.nc --clobber
    # Tracker program gettrk_main.f v3.9a can't handle variables with a vertical dimension.
    # It assumes a variable has only one pressure level.
    # delete unneeded variables.
    ncks -O -C -x -v z_isobaric,z_iso_levels,temperaure_200hPa tmp.nc tmp2.nc
    mv tmp2.nc tmp.nc

    # Instead of using meanT_500_300 from MPAS diagnostics just average the 5 levels with ncwa.
    # That is actually what hwrf_tave.exe does. (none of the dp-weighting used in MPAS diagnostics).
    # Extract 3-D temperature and 1-D pressure coordinate variable - makes ncrename and ncwa run faster
    ncks -v t_isobaric,t_iso_levels tmp.nc t_isobaric.nc
    # Use ncwa to average between 30000 Pa and 50000 Pa. Use decimal point so it knows they are real values not index values. 
    ncwa -a t_iso_levels -d t_iso_levels,30000.,50000. -O t_isobaric.nc tmean.nc
    # Change name of t_isobaric variable to tmean
    ncrename -v t_isobaric,tmean -O tmean.nc
    ncks -A -v tmean tmean.nc tmp.nc
    rm t_isobaric.nc tmean.nc

    # Define vertical level type of MSLP and convert from hPa to Pa (x100)
    # To see if it is needed, ncdump variable 'mslp' and look for numbers starting 
    # with 10 or 9 followed by 2 digits, then a decimal,
    # then another digit.  This matches pressure in hPa in mslp dump. Pressure
    # in Pa should not match. 2013 diagnostics mslp output was in hPa; 2014 in Pa.
    ncdump -v mslp ../diag*.${initdate}_${h}.00.00.nc | tail | grep -P " (10|9)\d\d\.\d"
    if ($status == 0) then
        ncap2 -O -s 'mslp=100.*mslp;' all.nc mslp.nc
        ncks -A mslp.nc tmp.nc
        rm mslp.nc
    endif

    # Set reference time of relative time axis and set base units to hours.
    cdo -O -setreftime,$initdate,${h}:00:00,1hour tmp.nc tmp2.nc

    # Bypass longitude problem by outputting 0-360 with mpas_to_latlon. - Jan 2019
    # gettrk handles flipped latitude now. 
    # This is a mystery. If I shift the order of the variables in memory from 
    # -180-180 to 0-360, I don't need to flip latitude. 
    # I fiddled with this to get atcf output with 'W' and 'E' values, not just 'E'. 
    # Before it was outputting just 'E' values and negative values. Not good. 
    # TODO: output latlon files with 0-360 longitude from mpas_to_latlon (avoid this step)
    #ncks -O --msa -d lon,0.,180. -d lon,-180.,-0.25 tmp2.nc tmp.nc # TODO: avoid hard-coded -0.25
    #ncap2 -O -s 'where(lon < 0) lon=lon+360' tmp.nc tmp2.nc
    # Flip latitude dimension.
    # cdo invertlat corrupts t_iso_levels for some reason
    # cdo -O invertlat tmp2.nc tmp.nc
    # Don't need to flip latitude dimension. gettrk seems to handle it. min and max lat are okay, etc...
    # actually I did need to flip latitude for v3.9a. 
    #ncpdq -O -a -lat tmp2.nc tmp.nc

    mv tmp2.nc $out
    rm tmp.nc
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
    echo removing $lock from $workdir
    rm $lock
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
  trkrinfo%northbd=50,
  trkrinfo%southbd=-5,
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
        rundescr='$mp',
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
verb=2
/
NL

printf "Running tracker..."
$BINDIR/gettrk.exe < namelist > log
echo done
if ($status != 0) then
    tail log
    exit 1
endif

cp fort.15 fort.61 fort.62 fort.64 fort.68 fort.69 namelist log $tcvitalsfile $trackertype
if($trackertype == 'tcgen') mv fort.66 fort.67 fort.77 $trackertype
echo removing lock file $lock
rm $lock

# jump here if --justplot option 
PLOT:
# Between tracker and tcgen, we prefer tcgen online.
set atcf_with_warmcore_column=fort.64
set toserver=""
if($trackertype == 'tcgen') then
    set atcf_with_warmcore_column=fort.66
    set toserver="--to_server"
endif

# Read tracks and plot with script
# -f force overwrite 
# plot_atcf.py --to_server rsyncs to web server
# The 1st basin sets the lat/lon grid tick interval for all basins.
if ($user != ahijevyc) echo only user ahijevyc can rsync to server
set echo
# TODO: allow cp basin. when I added cp to the basin list, it collapsed the global plot longitude range to zero.
python ~ahijevyc/bin/plot_atcf.py $workdir/gfdl_tracker/$trackertype/$atcf_with_warmcore_column global al wp ep io --stormname allstorms \
    --project $mp --origmesh --diagdir $idir/$mp/$yyyymmddhh --initfile $idir/$mp/$yyyymmddhh/init.nc --force_new $toserver

