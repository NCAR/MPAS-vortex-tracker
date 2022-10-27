#!/bin/csh
#
#BSUB -P P64000080
#BSUB -n 1
#BSUB -q geyser
#BSUB -W 1:45
##BSUB -w "done(run_mpas_to_latlon.csh)"
#BSUB -e /glade/scratch/ahijevyc/mpas_ll_GRIB1.%J.err
#BSUB -o /glade/scratch/ahijevyc/mpas_ll_GRIB1.%J.out
#BSUB -u ahijevyc


# convert hwt2017 (hourly) MPAS lat/lon files to GRIB for GFDL vortex tracker
# remember to interpolate to lat/lon with run_mpas_to_latlon.csh first
# For example:
#   ~ahijevyc/bin/run_mpas_to_latlon.csh 2017090700/ens_1 -w hwt2017 -t /glade/work/ahijevyc/tracking_gfdl/hwt2017_fields_to_interpolate
#
# usage
# mpas_hwt2017_ll_GRIB.csh -w hwt2017 [-i /glade/scratch/ahijevyc] yyyymmddhh

module load cdo
module load nco
module load python # for forecast_hour_links.py

if ( $HOSTNAME =~ cheyenne* ) then
    set BINDIR = /glade/u/home/ahijevyc/bin_cheyenne
else
    set BINDIR = /glade/u/home/ahijevyc/bin
endif 

set ymd_in=`date -u +%Y%m%d`
set mp=mpas\*
set trackertype=tcgen # tracker is traditional and only tracks stuff mentioned at t=0.  tcgen does genesis too
set idir=/glade/scratch/mpasrt
set debug=0

while ("$1" != "")
	if ("$1" =~ 20[0-9][0-9][01][0-9][0-3][0-9]*) set ymd_in="$1"
	if ("$1" == "-i") then # optional argument -i determines parent working directory. 
		shift
		set idir="$1"
	endif
	if ("$1" == "-w") then # optional argument -w determines working directory. added this because the # of potential working directory patterns is getting out of hand.
		shift
		set mp="$1"
	endif
	if ("$1" == "-t") then # optional argument -t can determine trackertype. (tracker or tcgen)
		shift
		set trackertype="$1"
	endif
	if ("$1" == "-d") set debug=1
	shift
end

if ("$1" != "") then
	echo unknown argument $1
	exit
endif

set tmpnc=mpas_ll_GRIB1.tmp.nc

set dxdetails=_0.500deg_025km
umask 2
if ($debug) set echo

foreach d ($idir/$mp/$ymd_in*)
	cd $d/latlon$dxdetails
	set refdate=`echo $d|sed -e 's,.*/\([12][09][0-9][0-9][01][0-9][0123][0-9][012][0-9]\).*,\1,'`
	set bcc=`echo $refdate|cut -c1-2`
	set byy=`echo $refdate|cut -c3-4`
	set m=`echo $refdate|cut -c5-6`
	set dd=`echo $refdate|cut -c7-8`
	set h=`echo $refdate|cut -c9-10`
	set ymd=$bcc$byy$m$dd
	set refdate=$bcc$byy-$m-$dd
	mkdir -p gfdl_tracker
	# used to use simple wildcard, but it matched hours between multiples of 3
	# Hourly files were present for random dates, like mpas3/2013092700.
	set if=0
	set fort15=gfdl_tracker/fort.15
	if (-e $fort15) rm $fort15
	foreach f (`ls diag*.$bcc??-??-??_??.00.00${dxdetails}.nc|sort` )

		# Append line to fort.15 (index of file starting with "0001" and forecast lead time in minutes) 
		# Used to be at end of foreach block but it needs to occur with every iteration of the diagnostics
		# file loop, before any "continue" clause.
		# Use forecast_hour_links.py to output the number of minutes. Allows for missing diag files. 
		set fmin=`~ahijevyc/bin/forecast_hour_links.py -m $f`
		printf '%04d %05d\n' `expr $if + 1` $fmin >> $fort15
		@ if++


		set out=gfdl_tracker/diag.$mp.$ymd$h.f`printf '%05d' $fmin`
		if (-s $out) then
			continue
		endif
		# Used to include RELV, but GFDL vortex tracker derives RELV from ABSV.
        # If you provide RELV directly, it doesn't know what to do with it. 
		# It will derive RELV from u/v. 
		ncap2 -O -s 'defdim("lv_ISBL0",3);lv_ISBL0[$lv_ISBL0]={85000,70000,50000};\
			u[time,$lv_ISBL0,lat,lon]=0.;\
			v[time,$lv_ISBL0,lat,lon]=0.;\
			u(:,0,:,:)=uzonal_850hPa;v(:,0,:,:)=umeridional_850hPa;\
			u(:,1,:,:)=uzonal_700hPa;v(:,1,:,:)=umeridional_700hPa;\
			u(:,2,:,:)=uzonal_500hPa;v(:,2,:,:)=umeridional_500hPa;\
			geopotential_height[time,$lv_ISBL0,lat,lon]=0.;\
			geopotential_height(:,0,:,:)=height_850hPa;\
			geopotential_height(:,1,:,:)=height_700hPa;\
			' $f $tmpnc

        if ($status != 0) exit
        # Extract netCDF fields Z, U, V. Ignore others. They will be added later.
        # ncks extracts and writes in alphabetical order by default. I listed these in 
        # alphabetical order to be consistent, but changing the order on the command line won't change the
        # output order (unless -a option is used).
        # Order is important for cdo, since cdo uses ordinal numbers -1, -2, etc. to refer to variables, not names.
		ncks -O -v geopotential_height,u,v $tmpnc $tmpnc
        # Convert to grib.
        cdo -f grb chparam,-1,7.2,-2,33.2,-3,34.2 $tmpnc zuv.grb

        touch tave.grb 

		# Define vertical level type of MSLP and convert from hPa to Pa (x100)
		ncks -O -v mslp $f $tmpnc
		# look for numbers starting with 10 or 9 followed by 2 digits, then a decimal,
		# then another digit.  This matches pressure in hPa in mslp dump. Pressure
		# in Pa should not match. 2013 diagnostics mslp output was in hPa; 2014 in Pa.
		ncdump -v mslp $tmpnc | tail | grep -P " (10|9)\d\d\.\d"
		if ($status == 0) ncap2 -O -s 'mslp=100.*mslp;' $tmpnc $tmpnc

        #cdo -f grb chparam,-1,2.2 -setltype,102 $tmpnc u.grb
        # With hwrf_gettrk v3.8a mslp must be 130 (membrane Eta reduction SLP), not 2. 
        cdo -f grb chparam,-1,130.2 -setltype,102 $tmpnc u.grb

		# Define grib level type of u10 and v10
		ncks -O -v u10,v10 $f $tmpnc
		cdo -f grb chparam,-1,33.2,-2,34.2 -setltype,105 -setlevel,10 $tmpnc v.grb

		# Merge grib records.(used to skip tave.grb for 2013).
		cat u.grb zuv.grb v.grb tave.grb > all.grb

		# Flip latitude dimension. Set reference time of relative time axis and set base units to hours.
		cdo -O invertlat -setreftime,$refdate,${h}:00:00,1hour all.grb $out

        # grb index file for tracker
		$BINDIR/grbindex.exe $out $out.ix

		# Clean up
		rm $tmpnc zuv.grb u.grb v.grb tave.grb all.grb
		

	end

	# Enter directory for GFDL vortex tracker files
	cd gfdl_tracker

	# Download combined TC vitals for requested year, save in fort.12.
	# Extract the date and hour that match the current init time
	# Consider using /glade/scratch/wrfrt/mpas/tracks/message*.$ymd instead. I found these Nov 2014 but don't know when they come in.
	# Consider using --fiorino option. 
	~ahijevyc/bin/wget_tcvitals.csh $ymd $h fort.12
	if($trackertype == 'tracker' && ! -s fort.12 )then
		echo "No storms at initial time. Exiting $0."
		touch no_tracker_storms
		exit
	endif  

	# Make empty fort.14 file.
	touch fort.14

	set phs=both
	cat <<NL > namelist

&datein
  inp%bcc=$bcc,
  inp%byy=$byy,
  inp%bmm=$m,
  inp%bdd=$dd,
  inp%bhh=$h,
  inp%model=1,
  inp%lt_units='hours',
  inp%file_seq='multi',
  inp%modtyp='global',
  inp%nesttyp='fixed'
/
&atcfinfo
  atcfnum=0,
  atcfname='MPAS',
  atcfymdh=$bcc$byy$m${dd}$h,
  atcffreq=100
/
&trackerinfo
  trkrinfo%westbd=-180,
  trkrinfo%eastbd=180.,
  trkrinfo%northbd=55,
  trkrinfo%southbd=0,
  trkrinfo%type='$trackertype',
  trkrinfo%mslpthresh=0.0015,
  trkrinfo%v850thresh=1.5000,
  trkrinfo%gridtype='global',
  trkrinfo%contint=100.0,
  trkrinfo%out_vit='y'
/
&phaseinfo
  phaseflag='n',
  phasescheme='$phs'
  wcore_depth=1.0
/
&structinfo
  structflag='y',
  ikeflag='n'
/
&fnameinfo
  gmodname='diag',
  rundescr='$mp',
  atcfdescr=''
/
&waitinfo
  use_waitfor='n'
/
&verbose
   verb=2
/
NL

	# Run tracker. Create tracks.
    $BINDIR/hwrf_gettrk.exe < namelist > log
	mkdir -p $trackertype
	cp fort.1[25] fort.61 fort.62 fort.64 fort.68 fort.69 fort.7? namelist log $trackertype
	if($trackertype == 'tcgen') mv fort.66 fort.67 $trackertype

end
