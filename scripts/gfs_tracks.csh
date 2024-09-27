#!/bin/csh
#
# Find GFS forecasts and run vortex tracker on them.
#
# For ACC, use p/work/ahijevyc/ncl/ACC.multi.ncl
#
# with no argument, assume the current date
# Otherwise the 1st argument can be an old date
#
# USAGE: gfs_tracks.csh [yyyymmddhh] [-t trackertype] [-gs grid_spacing]
#
#

# If you are looking to do ACC for 2016 and later, use run_fcst-init.csh
# and loop through grib2 files in /glade/collections/rda/data/ds084.1/

source ~ahijevyc/.tcshrc
module load grads # vertical average temperature (Tave) in grib2 file.
module load cdo # flip latitude dimension of grads output (Tave)
set BINDIR=~ahijevyc/bin
# to get wgrib2 on cheyenne, module load grib-bins
uname -a
module load grib-bins
if ($NCAR_HOST == cheyenne) then
    set BINDIR=~ahijevyc/bin_cheyenne
endif

# Default is today 00 UTC
set ymdh=`date -u +%Y%m%d`00
# set do_tracker to 1 to do vortex tracker by default.
# set do_tracker to 0 if you just want to download grib2 files. 
set do_tracker=1
set trackertype=tcgen
set gs=0.25

# Loop through arguments
while ("$1" != "")
	# added possibility of ymdh=20150930/ec
	if ("$1" =~ 20??????[01][02]*) set ymdh="$1" # added possibility of 12UTC 
	if ("$1" == "-nt") set do_tracker=0 # "nt" stands for "no tracker"
	if ("$1" == "-t") then # optional argument -t can determine trackertype. (tracker or tcgen)
		shift
		set trackertype="$1"
	endif
	if ("$1" == "-gs") then # grid spacing (0.5, 0.25)
		shift
		set gs="$1"
	endif
	shift
end

# Check for leftover arguments
if ("$1" != "") then
	echo unknown argument $1
	exit
endif


set year=`echo $ymdh|cut -c1-4`
set mm=`echo $ymdh|cut -c5-6`
set h=`echo $ymdh|cut -c9-10`
set ymd=`echo $ymdh|cut -c1-8`

set dir=/glade/scratch/ahijevyc/GFS/$ymdh/gfdl_tracker
mkdir -p $dir
cd $dir

set lock=.lock
if (-e $lock) then
	echo found lock file $dir/$lock
	exit
endif
touch $lock
# 3-digit forecast hours
# switched from 3 to 6-hour interval Nov 20 2014.
set dt=6
if ($do_tracker == 0) set dt=24
if ($do_tracker == 0 && $h == "12") set dt=12
set ff_list="`seq -w 000 $dt 192`"

foreach ff ($ff_list)
	set grb2=gfs.t${h}z.pgrb2.0p50.f${ff}

	set expected_records=358
	if ($year == 2014) then
		set expected_records=358
		if ($ff == 192) set expected_records=360 # analysis has fewer records than forecast
	endif
	if ( $ff == 000 ) set expected_records=315 # analysis has fewer records than forecast

	# on 0.25 deg grid
	if ( $gs == 0.25 ) then
		set grb2=/glade/collections/rda/data/ds084.1/$year/${ymd}/gfs.0p25.${ymdh}.f${ff}.grib2
		set expected_records=366
		if ($ff == 000) set expected_records=322 # analysis has fewer records than forecast
		if ($year >= 2016) then
			set expected_records=415
			if ($ff == 000) set expected_records=352 # analysis has fewer records than forecast
		endif
	endif

    # 2 records added to grib2 file after 20170719.
    if ( $ymd > 20170719 ) then
        set expected_records=`expr $expected_records + 2`
    endif

	echo expect $expected_records GRIB records from fhr $ff

	if ( ! -s $grb2 ) then
		echo $grb2 not found. Looking elsewhere.
		if ($ymdh >= `date +%Y%m%d`$h) then
		    echo Get GFS grib2 file from NCEP server. Only works with recent dates. 
			set grb2=gfs.t${h}z.pgrb2.0p25.f${ff}
			wget -nc http://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.${ymdh}/$grb2 --no-verbose
		else
			if ( $gs == 0.25 ) then
				echo "0.250° $grb2 not found. Skipping to next fcst time."
				continue
			endif
			# The NCDC archive can be used for 0.500deg and older dates.
			set grb2=gfs_4_${ymd}_${h}00_$ff.grb2 # grab individual hour
			wget -c -nc "ftp://nomads.ncdc.noaa.gov/GFS/Grid4/$year$mm/$ymd/$grb2" --no-verbose
		endif
	endif
	echo check integrity of AND number of records in GRIB2 file
	set nr=`wgrib2 $grb2 | wc -l`
	echo found $nr records
	if ( $nr != $expected_records ) then
		echo expected $expected_records
		exit
		if (-e $grb2) rm -v $grb2 # Sadly, there are corrupt grib2 files in ncdc and rda achive
		echo "try hsi:/RAPDMG/grib/GFS004/$ymd/${ymd}_i${h}_fxx_GFS004.grb2.gz on HPSS. see /glade/scratch/tomjr/gfs-isaac/getgfs_hpss.sc"
		echo skipping to next fcst time
		continue
	endif

	set fminutes=`echo $ff \* 60 | bc`
	set grb1=gfs.t${h}z.$ymdh.f`printf '%05d' $fminutes`
    
	if ( $do_tracker == 1 ) then
        if (! -w $grb2) then
            # If file is not writable copy to local directory.
            # You need a writable file so we may append Tave to it.
            # Just bring over variables needed for tracking
		    wgrib2 $grb2 -s |egrep '(:TMP:|:UGRD:|:VGRD:|:HGT:|:MSLET|:LAND)' | egrep -v '(:PV=|m above mean|trop|max wind|sigma|:30-0 mb)' | wgrib2 -i $grb2 -grib $grb1 
        endif
		# Make sure GRIB file has correct lat grid spacing. It may be an old one.
		set thisgs=`wgrib2 -V $grb1|head|grep lat|grep " by "|awk '{print $6}'`
		if (`echo "$gs != $thisgs"|bc`) echo requested $gs but got $thisgs
		# Make sure GRIB file has correct lon grid spacing. It may be an old one.
		set thisgs=`wgrib2 -V $grb1|head|grep lon|grep " by "|awk '{print $6}'`
		if (`echo "$gs != $thisgs"|bc`) echo requested $gs but got $thisgs
		set nr=`wgrib2 $grb1 | wc -l`
		# could be one more with TAVE added
		echo found $nr GRIB records. 

		# Add Tave, if necessary.
        # As of Nov 20, 2018, I haven't got this to work. cdo may not handle GRIB2.
		# Is 300-500 hPa vertically averaged temperature already in $grb1?
		wgrib2 $grb1 | grep TMP | grep "401 mb" > /dev/null
		if ($status != 0) then
            if (0) then
                # Calculate Tave with tave.exe.
                # hwrf_tave.exe made empty GRIB2 file. emailed Biswas Aug 31, 2017 for help.
                wgrib2 -v2 $grb1 | grep ":TMP" | grep ":lvl1=(100,"| wgrib2 -i $grb1 -grib fort.11
                # tave.exe expects input fort.11 and index file fort.31
                # Prepare index file for tavg.exe.
                $BINDIR/grbindex.exe fort.11 fort.31
                # fort.16 needed by tave.exe (Text file containing number of input pressure levels)
                wgrib2 fort.11 | wc -l > fort.16
                # vertically average temperature in 500-300 hPa layer
                echo "&timein ifcsthour=$ff,iparm=11,gribver=2/" | $BINDIR/tave.exe #> /dev/null
                if ($status != 0 ) then
                    echo problem with tave.exe
                    goto unlock
                endif
                # merge average temperature in 500-300 hPa layer and grib1 file
                cat fort.51 >> $grb1
                echo added tave to $grb1
                rm fort.51 fort.11
            else
                # Until hwrf_tave.exe is fixed and can write GRIB2...
                # Hack 300-500mb Tave with GRADS
                if ($gs != 0.25) then
                    echo "GrADS assumes 0.25deg grid spacing to calculate Tave. Exiting."
                    exit
                endif
                if (! -s grb.template) wgrib2 $grb1 -d 1 -lola 0:1440:0.25 -90:721:0.25 grb.template grib
                g2ctl $grb1 > $grb1.ctl
                gribmap -i $grb1.ctl
                # Use GrADS to average in the vertical and shave off repeated dateline longitude.
                grads -blcx "$HOME/bin/tave.gs $grb1.ctl $ymdh $ff"
                #cdo invertlat tave.grb2 t.grb2
                # cdo doesn't recognize grib2? Use wgrib2 instead.
                wgrib2 tave.grb -new_grid_winds earth -new_grid latlon 0:1440:0.25 90:721:-0.25 t.grb
                cat t.grb >> $grb1
                rm t.grb 
            endif
		endif
		# Prepare index file for GFDL vortex tracker.
		$BINDIR/grb2index.exe $grb1 $grb1.ix

	endif

end # ff_list foreach loop

if ($do_tracker == 1 && -s $grb1.ix) then # make sure at least one grib file exists. 

	# Prepare for GFDL vortex tracker.
	# Download combined TC vitals for requested year, save in fort.12.
	wget_tcvitals.csh $ymd $h fort.12
	if($trackertype == 'tracker' && ! -s fort.12) then
		echo "No storms to track at initial time $ymdh. Exiting $0." > no_tracker_storms
		goto unlock
	endif  

	# Make empty fort.14 file.
	touch fort.14

	# Write enumerated list of forecast times (in minutes) to fort.15
	# using all *.ix files as a basis.
	create_fort.15 *.ix


	# Extract date/time substrings for namelist
	set bcc=`echo $ymdh|cut -c1-2`
	set byy=`echo $ymdh|cut -c3-4`
	set bmm=`echo $ymdh|cut -c5-6`
	set bdd=`echo $ymdh|cut -c7-8`
	set bhh=`echo $ymdh|cut -c9-10`

	# Track TCs
cat <<NL > namelist
&datein
inp%bcc=$bcc,
inp%byy=$byy,
inp%bmm=$bmm,
inp%bdd=$bdd,
inp%bhh=$bhh,
inp%model=1,
inp%lt_units='hours',
inp%file_seq='multi',
inp%modtyp='global',
inp%nesttyp='fixed'
/
&atcfinfo
atcfnum=0,
atcfname='GFSO',
atcfymdh=$ymdh,
atcffreq=${dt}00
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
trkrinfo%out_vit='y',
trkrinfo%inp_data_type='grib',
trkrinfo%gribver=2
trkrinfo%g2_jpdtn=0,
trkrinfo%g2_mslp_parm_id=192,
/
&phaseinfo
phaseflag='y',
phasescheme='both'
wcore_depth=1.0
/
&structinfo
structflag='n',
ikeflag='n'
/
&fnameinfo
gmodname='gfs',
rundescr='t${h}z',
atcfdescr=''
/
&waitinfo
use_waitfor='n'
/
&verbose
verb=2
/

NL

    $BINDIR/gettrk.exe < namelist > log
    if ($status != 0) goto unlock

	# Save files to subdirectory (tracker, tcgen, or tcgen_0.50deg)
	set subdir=$trackertype
	# If 0.5-deg, prepend to subdirectory name.
	if($gs == 0.5)set subdir=${trackertype}_0.5deg
	mkdir -p $subdir
	cp fort.1[25] fort.61 fort.62 fort.64 fort.68 fort.69 fort.7? namelist log $subdir
	if($trackertype == 'tcgen') mv fort.66 fort.67 $subdir

	# Read tracks and plot with script
	# -f force overwrite 
    module load python
    ncar_pylib
    setenv PYTHONPATH /glade/u/home/ahijevyc/lib/python2.7
    # plot_atcf.py rsyncs to nova too.
	if($trackertype == 'tcgen') python ~ahijevyc/bin/plot_atcf.py $subdir/fort.66 --no-origmesh --force_new --to_server

	# Copy to Tom G.'s directory on MMM web server.
	#~ahijevyc/bin/Copy2WebTomjr.csh GFSO_${ymdh}_
endif

rm $lock
exit 0

unlock:
rm $lock
exit 1
