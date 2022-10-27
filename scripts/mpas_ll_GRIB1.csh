#!/bin/csh
#
# convert MPAS lat/lon files to GRIB1 for GFDL vortex tracker
# remember to interpolate to lat/lon with run_mpas_to_latlon.csh first
#
# Compliant with v3.9 of vortex tracker
#
# usage
# mpas_ll_to_GRIB1.csh [yyyymmddhh] [-w model] [-t (tracker|tcgen)] [-d] [-i workdir_parent] 

module load cdo
module load nco
module load python # for forecast_hour_links.py and ncar_pylib

set BINDIR = /glade/u/home/ahijevyc/bin
printenv NCAR_HOST
if ( $NCAR_HOST =~ cheyenne* ) then
    set BINDIR = /glade/work/ahijevyc/standalone_gfdl-vortextracker_v3.9a/trk_exec
    set BINDIR = /glade/u/home/ahijevyc/bin_cheyenne
endif 
if ( $NCAR_HOST =~ dav ) then
    set BINDIR = /glade/u/home/ahijevyc/bin_dav
endif 

set ymd_in=`date -u +%Y%m%d`
set mp=mpas\*
set trackertype=tcgen # tracker is traditional and only tracks stuff mentioned at t=0.  tcgen does genesis too
set idir=/glade/scratch/mpasrt
set debug=0

while ("$1" != "")
    # Loop thru multiple dates by specifying a date prefix like 2018060 instead of 2018060100
	if ("$1" =~ 20[0-9][0-9][01][0-9][0-3]*) set ymd_in="$1"
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

if ($debug) set echo

if ("$1" != "") then
	echo unknown argument $1
	exit
endif

set tmpnc=mpas_ll_GRIB1.tmp.nc

set dxdetails=_0.500deg_025km
umask 2

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
	mkdir -p gfdl_tracker/$trackertype
	# used to use simple wildcard, but it matched hours between multiples of 3
	# Hourly files were present for random dates, like mpas3/2013092700.
	set if=0
	set fort15=gfdl_tracker/fort.15
	if (-e $fort15) rm $fort15
	foreach f (`ls diag*.$bcc??-??-??_0[06].00.00${dxdetails}.nc \
	           diag*.$bcc??-??-??_1[28].00.00${dxdetails}.nc|sort` )

		# Append line to fort.15 (index of file starting with "0001" and forecast lead time in minutes) 
		# Used to be at end of foreach block but it needs to occur with every iteration of the diagnostics
		# file loop, before any "continue" clause.
		# Use forecast_hour_links.py to output the number of minutes. Allows missing diag files. 
		set fmin=`~ahijevyc/bin/forecast_hour_links.py -m $f`
		printf '%04d %05d\n' `expr $if + 1` $fmin >> $fort15
		@ if++


		set out=gfdl_tracker/diag.$mp.$ymd$h.f`printf '%05d' $fmin`
		if (! -s $out) then
            # Used to include RELV, but GFDL vortex tracker derives RELV from ABSV.
            # If you provide RELV directly, it doesn't know what to do with it. 
            # It will derive RELV from u/v. 
            ncap2 -O -s 'defdim("pressure",3);pressure[$pressure]={85000,70000,50000};\
                u[time,$pressure,lat,lon]=0.;\
                v[time,$pressure,lat,lon]=0.;\
                u(:,0,:,:)=uzonal_850hPa;v(:,0,:,:)=umeridional_850hPa;\
                u(:,1,:,:)=uzonal_700hPa;v(:,1,:,:)=umeridional_700hPa;\
                u(:,2,:,:)=uzonal_500hPa;v(:,2,:,:)=umeridional_500hPa;\
                ' $f $tmpnc
            # rename T and Z pressure levels variables something cdo understands with GRIB
            ncrename -O -d nIsoLevelsT,lv_ISBL0 -d nIsoLevelsZ,lv_ISBL1 -v t_iso_levels,lv_ISBL0 -v z_iso_levels,lv_ISBL1 $tmpnc

            # Extract netCDF fields T, U, V, and Z. Ignore others. They will be added later.
            # ncks extracts and writes in alphabetical order by default. I listed these in 
            # alphabetical order to be consistent, but changing the order on the command line won't change the
            # output order (unless -a option is used).
            # Order is important for cdo, since cdo uses ordinal numbers -1, -2, etc. to refer to variables, not names.
            ncks -O -v t_isobaric,u,v,z_isobaric $tmpnc $tmpnc
            # Convert to grib.
            cdo -f grb chparam,-1,11.2,-2,33.2,-3,34.2,-4,7.2 $tmpnc tuvz.grb

            # 200 mb height
            ncks -O -v height_200hPa $f $tmpnc
            cdo -f grb chparam,-1,7.2 -setltype,100 -setlevel,200 $tmpnc z200.grb

            # Instead of using meanT_500_300 from MPAS diagnostics just average the 5 levels with ncwa.
            # That is actually what hwrf_tave.exe does. (none of the dp-weighting used in MPAS diagnostics).  
            ncwa -a nIsoLevelsT -v t_isobaric -O $f $tmpnc
            cdo -f grb chparam,-1,11.2 -setltype,100 -setlevel,401 $tmpnc tave.grb

            # Define vertical level type of MSLP and convert from hPa to Pa (x100)
            ncks -O -v mslp $f $tmpnc
            # To see if it is needed, ncdump variable 'mslp' and look for numbers starting 
            # with 10 or 9 followed by 2 digits, then a decimal,
            # then another digit.  This matches pressure in hPa in mslp dump. Pressure
            # in Pa should not match. 2013 diagnostics mslp output was in hPa; 2014 in Pa.
            ncdump -v mslp $tmpnc | tail | grep -P " (10|9)\d\d\.\d"
            if ($status == 0) ncap2 -O -s 'mslp=100.*mslp;' $tmpnc $tmpnc

            #cdo -f grb chparam,-1,2.2 -setltype,102 $tmpnc u.grb
            # With hwrf_gettrk v3.8a mslp must be 130.2 (membrane Eta reduction SLP), not 2.2. 
            cdo -f grb chparam,-1,130.2 -setltype,102 $tmpnc u.grb

            # Define vertical level type and vertical level of u10 and v10 in a standard
            # way recognized by grib readers.
            ncks -O -v u10,v10 $f $tmpnc
            cdo -f grb chparam,-1,33.2,-2,34.2 -setltype,105 -setlevel,10 $tmpnc v.grb

            # Merge grib records.(used to skip tave.grb for 2013).
            cat tave.grb tuvz.grb u.grb v.grb z200.grb > all.grb

            # Flip latitude dimension. Set reference time of relative time axis and set base units to hours.
            cdo -O invertlat -setreftime,$refdate,${h}:00:00,1hour all.grb $out
            # Clean up
            rm $tmpnc tave.grb tuvz.grb u.grb v.grb z200.grb all.grb
            
		endif

        # grb index file for tracker
		$BINDIR/grbindex.exe $out $out.ix


	end

	# Enter directory for vortex tracker files
	cd gfdl_tracker

    # put lock file in working directory of tracker. 
    # Don't want more than 1 tracker running in the same directory.
    set lock=.mpas_ll_gettrk.lock
    if (-e $lock) then
        echo found lock file $lock in $d 
        continue
    endif
    touch $lock

	# Download combined TC vitals for requested year, save in 
    # tcvit_rsmc_storms.txt (tracker, formerly fort.12).
    # or tcvit_genesis_storms.txt (tcgen, formerly fort.12)- couldn't get this to work with tcvit_rsmc_storms.txt existing.
    set tcvitalsfile=tcvit_rsmc_storms.txt
    #if ($trackertype == 'tcgen') set tcvitalsfile=tcvit_genesis_storms.txt
	# Extract the date and hour that match the current init time
	# Consider using --fiorino option. 
	~ahijevyc/bin/wget_tcvitals.csh $ymd $h $tcvitalsfile
	if($trackertype == 'tracker' && ! -s $tcvitalsfile )then
		echo "No storms at initial time $ymd$h. Exiting $0." > no_tracker_storms
        cp no_tracker_storms $trackertype
        echo removing $lock from $d
		rm $lock
        continue
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
  atcffreq=300
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
  trkrinfo%inp_data_type='grib',
  trkrinfo%gribver=1,
  trkrinfo%g1_mslp_parm_id=130,
  trkrinfo%g1_sfcwind_lev_typ=105,
  trkrinfo%g1_sfcwind_lev_val=10
/
&phaseinfo
  phaseflag='y',
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
&netcdflist
/
&parmpreflist
/
&verbose
 verb=3
/
NL

	# Run tracker. Create tracks.
    $BINDIR/gettrk.exe < namelist > log
    if ($status != 0) then
        tail log
        exit 1
    endif
	cp fort.15 fort.61 fort.62 fort.64 fort.68 fort.69 fort.7? namelist log $tcvitalsfile $trackertype
	if($trackertype == 'tcgen') mv fort.66 fort.67 $trackertype

	# Read tracks and plot with script
	# -f force overwrite 
    module load python
    ncar_pylib
    setenv PYTHONPATH /glade/u/home/ahijevyc/lib/python2.7
    # plot_atcf.py rsyncs to nova too.
	if($trackertype == 'tcgen') python ~ahijevyc/bin/plot_atcf.py $trackertype/fort.66 --force_new --to_server

    echo removing $lock from $d
    rm $lock
end
