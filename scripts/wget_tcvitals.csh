#!/bin/csh
umask 2
if ("$1" == "") then
	echo "USAGE: $0 [-fiorino] yyyymmdd hh outfile"
	echo "(outfile used to be fort.12)"
	exit
endif

if ("$1" == "--fiorino" || "$1" == "-fiorino") then
	# From Mike Fiorino email (Mar 4, 2016)
	# ...regenerated 'full' tcvitals (includes invests or 9X or pre/potential TCs or pTCs)
	# ...these tcvitals are more complete than any available from operations.
	# I take special and greater care of the pTC a- and b-decks than either JTWC or NHC.
	#see: http://ruc.noaa.gov/fiorino/tc/tcvitals-2014.tgz or http://ruc.noaa.gov/fiorino/tcvitals/tcvitals.YYYYMMDDHH.txt
	# where YYYYMMDDHH is the 'date-time-group' or year,month,day,hour
	set fiorino
	shift
endif


set ymd=$1
if ($ymd !~ [12]???????) then
	echo "bad ymd: '$ymd'"
	exit
endif
if ("$3" == "")  then
	echo "3rd argument must be output file name"
	exit
endif

# 2 digit hour
set h=`printf %02d $2`
# output filename 
set out=$3
# Download combined TC vitals for requested year
set yyyy=`echo $ymd|cut -c1-4`
set tmpfile=$TMPDIR/wget_tcvitals.$$

if ($?fiorino) then
	wget https://ruc.noaa.gov/fiorino/tc/tcvitals/tcvitals.$ymd$h.txt --no-verbose
	mv -v tcvitals.$ymd$h.txt $out
else
	set tcv=combined_tcvitals.$yyyy.dat
	# -N retrieve new file if local file is older
	wget http://hurricanes.ral.ucar.edu/repository/data/tcvitals_open/$tcv --output-file=$tmpfile
	if ($status != 0) wget http://hurricanes.ral.ucar.edu/repository/data/tcvitals_open/$yyyy/$tcv --output-file $tmpfile
	# Extract the date and hour that match the current init time
	grep "$ymd ${h}00 " $tmpfile | tee $out
endif
if (-e $tmpfile) rm $tmpfile
touch $tmpfile # existing, but empty file

# In Jonathan's archive, sometimes the same storm is listed twice with different lat/lons. 
# For example, 2016090900 INVEST 94L
# This results in a t=0 track location of 0N, 0E from GFDL vortex tracker. 
# Not a problem with NHC's archive but Jonathan's is more timely. 
#
# Eliminate multiple instances of the same storm. Keep the last.
#  the last one is guaranteed by the "tail" command. If you wanted the
#  first, then use "head".
# Convert spaces to underscores
set uniqstorms=`cat $out | cut -c6-19 | tr " " _ | sort | uniq`
echo $uniqstorms
foreach f ($uniqstorms)
	# Convert underscores back to spaces
	set g=`echo $f | tr _ " "`
	# need double quotes on $g to protect spaces
	grep "$g" $out | tail -n 1 >> $tmpfile
end
mv -v $tmpfile $out

