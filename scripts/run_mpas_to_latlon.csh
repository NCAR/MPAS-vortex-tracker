#!/bin/csh
#
# Interpolate to lat lon grid and filter_radius_km.
# uses mpas_to_latlon, a f90 program adapted from Michael Duda.
# mpas_to_latlon was wired to derive height of model surfaces and vorticity but I commented it out. 
# now it just processes the fields in the standard input list.

# Used to break into mpas/mpas2 and 08/16 (4 parts) and run in parallel.
# But now that I added a lock file, just start the whole thing many times.
# The program will skip over completed times and partially-completed times.

# Used to find vortices, but commented this out. Now tracking is done by GFDL vortex tracker and
# handled by mpas_ll_GRIB1.csh. 

module load nco
umask 002

set debug=0
set delta=0.5
set lat0=-30
set lat1=60
set lon0=0
set filter_radius_km=25
set EXECUTABLE=~ahijevyc/bin/mpas_to_latlon

set current_date=`date -u +%Y%m%d`00
set ymdh=$current_date
set meshid=mpas_conv
set fields_to_interpolate=/glade/work/ahijevyc/tracking_gfdl/mpas_fields_to_interpolate.txt
set workdir=/glade/scratch/$user

while ("$1" != "")
    # Used to expect plain date directories.
	# Now you can append a subdirectory like ymdh=2015093000/ec or 2016092300/test
	if ("$1" =~ 20[0-9][0-9][01][0-9][0-3][0-9][01][02]*) set ymdh="$1" # added possibility of 12UTC 
	if ("$1" == "-d") then 
		set debug=1
	endif
	if ("$1" == "--delta") then # optional argument -d overrides default grid spacing 0.5deg
		shift
		set delta="$1"
	endif
	if ("$1" == "--executable") then 
		shift
		set EXECUTABLE="$1"
	endif
	if ("$1" == "--filter_radius_km") then # optional argument overrides '25' default
        # $filter_radius_km used to be the radius of the smoothing filter in grid points prior to July 2014. 
        # Now it is the radius in kilometers.
		shift
		set filter_radius_km="$1"
	endif
	if ("$1" == "-t") then # optional argument -t points to a text file with fields to interpolate
		shift
		set fields_to_interpolate="$1"
	endif
	if ("$1" == "--lat0") then # optional argument lat0 overrides default southern latitude
		shift
		set lat0="$1"
	endif
	if ("$1" == "--lat1") then # optional argument lat1 overrides default northern latitude
		shift
		set lat1="$1"
	endif
	if ("$1" == "--lon0") then # optional argument lon0 overrides default western longitude (gettrk likes 0)
		shift
		set lon0="$1"
	endif
	if ("$1" == "--mesh") then 
		shift
		set meshid="$1"
	endif
	if ("$1" == "-w") then # optional argument -w determines working directory.
		shift
		set workdir="$1"
	endif
	shift
end

if ("$1" != "") then
	echo unknown argument $1
	exit 1
endif

if (! -s $fields_to_interpolate) then
    echo "Couldn't find list of fields to interpolate: '$fields_to_interpolate'"
    exit 2
endif

set d=$workdir/$meshid/$ymdh
if (! -d $d) then
    echo "Did not find $d. Are working directory, mesh id, and initialization time correct?"
    echo "\tworking directory   $workdir"
    echo "\tmesh id             $meshid"
    echo "\tinitialization time $ymdh"
    exit 2
endif

cd $d

uname -n
if ( $debug ) then
    uname -a
    printenv
    set echo
endif

# check file existence. This file is dumped out after the last ncl job
# is finished.'
waitfile0:
set check_file=log.0000.out
# If running on the current date, wait until check_file is there. 
if ($ymdh =~ $current_date && ! -s $check_file) then
	echo in $d MPAS not finished. will try again in 10 minutes.
	sleep 600
	goto waitfile0
endif 

# Normally date directories have 3-hourly files, but some have hourly. 
# capture 6-hour multiples, but not 3-hour. Nov 20 2014
# just 6-h multiples Jun 2015
# For debugging sometimes Wei names them 'diag' instead of 'diagnostics'
# Don't zero pad. Let printf %02d do that. Otherwise it will consider it octal.
set dxdetails=`printf '_%5.3fdeg_%03dkm' ${delta} ${filter_radius_km}`
set odir=latlon$dxdetails
mkdir -p $odir/gfdl_tracker
if ($status != 0) then
    echo problem making output directory $odir/gfdl_tracker
    exit 2
endif
foreach f (diag*.20[123]?-??-??_0[06].00.00.nc diag*.20[123]?-??-??_1[28].00.00.nc )
    set ofile=$odir/$f
    # Only do this file if it doesn't exist already. or if $f is newer than $ofile.
    # `ls -t1 $f $ofile` lists files one on each line with newest first
    if (! -e $ofile || `ls -1t $f $ofile|head -n 1` == $f ) then
        set lock = $odir/.$f
        if (-e $lock) then
            echo found lock file $lock  moving on.
            continue
        endif
        touch $lock

        # limit vertical levels to these essential levels so mpas_to_latlon doesn't take so long.
        ncrename -d nIsoLevelsU,u_iso_levels -d nIsoLevelsT,t_iso_levels -d nIsoLevelsZ,z_iso_levels -O $f $lock
        ncks -d u_iso_levels,50000.,85000. -d t_iso_levels,30000.,50000. -d z_iso_levels,20000.,90000. -O $lock $lock

        set args = (`printf '%s %5.3f %02d %s %4.1f %4.1f %6.1f' $ofile $delta $filter_radius_km $meshid $lat0 $lat1 $lon0`)
        # Stored list of fields to interpolate in file. Oct 6, 2017. 
        echo "grep -v '#' $fields_to_interpolate | $EXECUTABLE $lock $args"
        grep -v '#' $fields_to_interpolate | $EXECUTABLE $lock $args

        if ($status == 2) then
            shift args # get rid of 1st element and move remaining down by 1
            echo "run echo precipw | $EXECUTABLE <path to init.nc> $odir/init.nc $args"
            exit
        endif

        if ($status != 0) then
            rm $ofile
            break
        endif

        # Append vertical coordinates to output file.
        # Don't append from original file $f because it will have all the vertical levels, not the filtered ones
        # filtered by the ncrename/ncks commands above. 
        foreach v (t z u)
            ncdump -h $f | grep ${v}_iso_levels >& /dev/null
            if ($status == 0) ncks -A -v ${v}_iso_levels $lock $ofile
        end


        rm $lock
    endif
end # diagnostics file

# moved a section to ~/bin/old/mpas-mpas2.csh in which I subtract fields in variable resolution (15-60km) MPAS from uniform resolution (15km) MPAS
# after interpolating to lat-lon grid.
