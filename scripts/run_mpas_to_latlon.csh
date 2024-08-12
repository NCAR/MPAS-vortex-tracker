#!/bin/csh
#
# Interpolate to lat lon grid and filter_radius_km.
# uses mpas_to_latlon, a f90 program adapted from Michael Duda.
# mpas_to_latlon was wired to derive height of model surfaces and vorticity but I commented it out. 
# now it just processes the fields in the standard input list.

# Used to find vortices, but commented this out. Now tracking is done by GFDL vortex tracker and
# handled by mpas_ll_GRIB1.csh. 

module load nco

set debug=0
set delta=0.5
set latmin=-5
set latmax=55
set lonmin=0
set lonmax=360
set filter_radius_km=25
set EXECUTABLE=$SCRATCH/MPAS-vortex-tracker/bin/mpas_to_latlon

set meshid=tk707_conus
set fields_to_interpolate=$SCRATCH/MPAS-vortex-tracker/scripts/mpas_fields_to_interpolate.txt
set idir=/glade/campaign/mmm/wmr/weiwang/cps/irma3/2020/tk707_conus
set workdir=$TMPDIR

while ("$1" != "")
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
	if ("$1" == "--latmin") then # optional argument latmin overrides default southern latitude
		shift
		set latmin="$1"
	endif
	if ("$1" == "--latmax") then # optional argument latmax overrides default northern latitude
		shift
		set latmax="$1"
	endif
	if ("$1" == "--lonmin") then # optional argument lonmin overrides default western longitude (gettrk likes 0)
		shift
		set lonmin="$1"
	endif
	if ("$1" == "--lonmax") then # optional argument lonmax overrides default maximum longitude
		shift
		set lonmax="$1"
	endif
	if ("$1" == "--mesh") then 
		shift
		set meshid="$1"
	endif
	if ("$1" == "-i") then # optional argument -i determines input directory.
		shift
		set idir="$1"
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

uname -n
if ( $debug ) then
    uname -a
    printenv
    set echo
endif


# Normally date directories have 3-hourly files, but some have hourly. 
# capture 6-hour multiples, but not 3-hour. Nov 20 2014
# just 6-h multiples Jun 2015
# For debugging sometimes Wei names them 'diag' instead of 'diagnostics'
# Don't zero pad. Let printf %02d do that. Otherwise it will consider it octal.
set dxdetails=`printf '_%5.3fdeg_%03dkm' ${delta} ${filter_radius_km}`
set odir=$workdir/$meshid/latlon$dxdetails
mkdir -p $odir/gfdl_tracker
if ($status != 0) then
    echo problem making output directory $odir/gfdl_tracker
    exit 2
endif
# used by mpas_ll_gettrk.csh
cat <<BOUNDS>bounds.csh
setenv latmin $latmin
setenv latmax $latmax
setenv lonmin $lonmin
setenv lonmax $lonmax
BOUNDS
foreach f ($idir/diag*.20??-??-??_0[06].00.00.nc $idir/diag*.20??-??-??_1[28].00.00.nc )
    set ofile=$odir/`basename $f`
    # Only do this file if it doesn't exist already. or if $f is newer than $ofile.
    # `ls -t1 $f $ofile` lists files one on each line with newest first
    if (! -e $ofile || `ls -1t $f $ofile|head -n 1` == $f ) then

        set args = (`printf '%s %5.3f %02d %s %4.1f %4.1f %6.1f %6.1f' $ofile $delta $filter_radius_km $meshid $latmin $latmax $lonmin $lonmax`)
        # Stored list of fields to interpolate in file. Oct 6, 2017. 
        echo "grep -v '#' $fields_to_interpolate | $EXECUTABLE $f $args"
        grep -v '#' $fields_to_interpolate | $EXECUTABLE $f $args

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
        foreach v (t z u)
            ncdump -h $f | grep ${v}_iso_levels >& /dev/null
            if ($status == 0) ncks -A -v ${v}_iso_levels $f $ofile
        end

    endif
end # diagnostics file

# moved a section to ~/bin/old/mpas-mpas2.csh in which I subtract fields in variable resolution (15-60km) MPAS from uniform resolution (15km) MPAS
# after interpolating to lat-lon grid.
