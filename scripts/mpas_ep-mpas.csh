#!/bin/csh

cd /glade/p/nmmm0024/mpas_ep

set dmpas="mpas_ep-mpas"
set outdir=/glade/work/ahijevyc/rmse_comparisons/$dmpas
mkdir -p $outdir
#goto dave

# A lot of this script is redundant with ~/bin/stddev.csh
# But stddev.csh deals with native MPAS mesh, not lat-lon grids.
#
# calls ~ahijevyc/bin/Hovmoller.csh file.nc lat0 lat1 

# 1) subtract mpas from mpas_ep
# Run all dates in parallel bsub jobs
# Remember to "escape" the dollar signs in the bsub input
# to delay interpretation.
#
foreach d ( 20*00/latlon_0.5*)
	if (`ls -1 $d/diag*.nc_minus_mpas|wc -l` == 41 || `dirname $d` >= 2014110400) continue
	echo submitting $d
	cat <<EOF | bsub
#!/bin/csh
#BSUB -P NMMM0024
#BSUB -n 1
#BSUB -q caldera
#BSUB -W 0:05
#BSUB -o /glade/scratch/ahijevyc/$dmpas.$$.out
#BSUB -J $dmpas
#BSUB -u ahijevyc 
foreach f ( $d/diag*.nc )
	set out=\${f}_minus_mpas
	if (! -s \$out && -s ../mpas/\$f) ncdiff \$f ../mpas/\$f \$out
end
EOF
sleep 1
end

# wait till all jobs are finished
wait_for _ep-mpas
echo finished subtracting mpas from mpas_ep


# link forecast hour filenames (fxxx.nc) to (mpas_ep - mpas) difference files  
foreach f (2014*00/latlon_0.5*)
	cd $f
	~ahijevyc/bin/forecast_hour_links.py *_minus_mpas
	if ($status != 0) echo problem making forecast hour links to mpas_ep-mpas in $f 
	cd ../..
end

# 2) get mean difference (mpas_ep - mpas) and average squared difference for each forecast hour.
# run all forecast lead times in parallel
foreach f (`seq -w 0 6 240`)
	echo submitting f$f
cat <<EOF | bsub
#!/bin/csh
#BSUB -P NMMM0024
#BSUB -n 1
#BSUB -q caldera
#BSUB -W 0:05
#BSUB -o /glade/scratch/ahijevyc/$dmpas.$$.out
#BSUB -J $dmpas
#BSUB -u ahijevyc 
foreach op (avg avgsqr)
set out=$outdir/$dmpas.\$op.f$f.nc
if (! -s \$out) then
nces -O --op_typ \$op 201*00/latlon_0.5*/f$f.nc \$out
# Append nfile variable
ncap2 -A -s "nfile[time]=`ls 201*00/latlon_0.5*/f$f.nc|wc -w`d;" \$out \$out
endif
end
EOF
end
wait_for _ep-mpas
echo finished averaging $dmpas difference and averaging the squared difference
# After all forecast hours are done, group and rename for mpas_basin_rmse.pro
# Aug 29-don't remember what I expected to do with this IDL program. . .it doesn't work on these nc files
# that have all the field variables together; it works on nc files with one field (u10_f-i.nc).
cd $outdir # added aug 29. not sure how it worked without it.
foreach op (avg avgsqr)
	ncrcat -O $dmpas.$op.f???.nc $dmpas.$op.nc
	# thought about removing individual forecast files, but messes up file checks earlier in script.
	# if ($status == 0) rm $dmpas.$op.f???.nc
end
ln -sf $dmpas.avg.nc $dmpas.ttl.nc


# Get normalization factor (time variance of uniform MPAS at each point).
# I used to say it was the standard deviation but I never take the square root.
cd /glade/p/nmmm0024/mpas
# in each latlon directory, link to diag* files
foreach f (20*00/latlon_0.5*)
	cd $f
	# Do not link to the "_minus_mpas" files. (require filename to end with ".nc")
	forecast_hour_links.py diag*.nc
	if ($status != 0) echo problem making forecast hour links in $f 
	cd ../..
end

# run all forecast lead times in parallel
# Used to include every forecast hour to smooth out noise from individual storms.
foreach f (`seq -w 0 6 240`)
	set avg=$outdir/mpas.avg.f$f.nc
	set avgsqr=$outdir/mpas.avgsqr.f$f.nc
	echo submitting f$f
cat <<EOF | bsub
#!/bin/csh
#BSUB -P NMMM0024
#BSUB -n 1
#BSUB -q caldera
#BSUB -W 0:05
#BSUB -o /glade/scratch/ahijevyc/$dmpas.$$.out
#BSUB -J $dmpas
#BSUB -u ahijevyc
# average for this forecast lead time 
if (! -s $avg) nces -O --op_typ avg 201*00/latlon_0.5*/f$f.nc $avg
# Deviation from average for each model run at this lead time
foreach f (2*00/latlon_0.5*/f$f.nc)
if (! -s \${f}_anom) ncbo --op_typ=- \$f $avg \${f}_anom
end
# Average of the squares
if (! -e $avgsqr) nces -y avgsqr 2*00/latlon_0.5*/f$f.nc_anom $avgsqr
EOF

end # forecast hour
# wait till all jobs are finished (search for substring "_ep-mpas". it is the only part that shows up in the bjobs output)
wait_for _ep-mpas
echo finished averaging the difference and averaging the squared difference 


ncrcat -O $outdir/mpas.avgsqr.f???.nc $outdir/mpas.avgsqr.nc
dave:
cd $outdir

# add u and v components of vector wind 
ncap2 -O -s "vector10=u10+v10;\
            vector_850hPa=uzonal_850hPa+umeridional_850hPa;\
            vector_700hPa=uzonal_700hPa+umeridional_700hPa;\
            vector_500hPa=uzonal_500hPa+umeridional_500hPa;" $dmpas.avgsqr.nc

$Hovm $dmpas.avgsqr.nc 0. 23.5
$Hovm $dmpas.avgsqr.nc 23.5 45.

#Here's some hard-wired stuff for 500hPa vector wind
set field=vector_500hPa
ncap2 -O -v -s "$field=sqrt($field/109.7299)" $dmpas.avgsqr.23.5-45..nc normalized_${field}_$dmpas.avgsqr.23.5-45..nc
ncap2 -O -v -s "$field=sqrt($field/33.53478)" $dmpas.avgsqr.0.-23.5.nc normalized_${field}_$dmpas.avgsqr.0.-23.5.nc
set field=height_200hPa
ncap2 -O -v -s "$field=sqrt($field/18893.291)" $dmpas.avgsqr.23.5-45..nc normalized_${field}_$dmpas.avgsqr.23.5-45..nc
ncap2 -O -v -s "$field=sqrt($field/692.52770)" $dmpas.avgsqr.0.-23.5.nc normalized_${field}_$dmpas.avgsqr.0.-23.5.nc



if (0) then 
	# Normalize global bias (get stddev and divide bias by it.)
	ncea -O -y sqrt $dmpas.avgsqr.nc $dmpas.stddev.nc
	ncbo -O --op_typ divide $dmpas.avg.nc $dmpas.stddev.nc normalized_$dmpas.avg.nc
	# Normalize global rmse  (divide avgsqr by variance and take the sqrt)
	ncbo -O --op_typ divide $dmpas.avgsqr.nc mpas.avgsqr.nc normalized_$dmpas.avgsqr.nc
	ncea -O -y sqrt normalized_$dmpas.avgsqr.nc normalized_$dmpas.avgsqr.nc

	# run a script called "Hovmoller.csh". It latitudinally-averages a lat-lon .nc file 
	# between lat0 and lat1 with cos(lat) weighting and outputs the new file's name.
	set Hovm=~ahijevyc/bin/Hovmoller.csh
	foreach lats ("0. 23.5" "23.5 45.")
		# First latitudinally-average the average of the squared differences in uniform mesh mpas.
		# This is the normalization factor.
		set variance=`$Hovm mpas.avgsqr.nc $lats`
		set stddev=`echo $variance|sed -e 's/avgsqr/stddev/'`
		echo variance=$variance stddev=$stddev
		ncea -O  -y sqrt $variance $stddev
		# Latitudinally average bias.
		set bias=`$Hovm $dmpas.avg.nc $lats`
		# Divide bias (mpas_ep-mpas.avg.nc) by the square root of variance.
		ncbo -O --op_typ divide $bias $stddev normalized_$bias
		# Latitudinally average mean squared difference.
		set msd=`$Hovm $dmpas.avgsqr.nc $lats`
		# Then divide by normalization factor.
		# Divide mean squared difference (mpas_ep-mpas.avgsqr.nc) by variance, then take square root.
		ncbo -O --op_typ divide $msd $variance normalized_$msd
		ncea -O -y sqrt normalized_$msd normalized_$msd
		echo normalized $lats
	end
	# Normalizing means the same range can be used for all fields in ncview.
endif

# Finally, to create Hovmollers, run $outdir/make_composite.csh

