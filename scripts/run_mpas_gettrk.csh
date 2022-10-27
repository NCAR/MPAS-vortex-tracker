#!/bin/csh

setenv TMPDIR /glade/scratch/$USER/temp
# Use unique id $$ to avoid conflict with other running dates.
set batchfile=$TMPDIR/mpas_gettrk.$$.sbatch

# 01:45:00 was too short after adding fancy background to plot_atcf.py
set walltime="02:00:00"


# set to 0 if you want to run on cheyenne share queue
# set to 1 if you want to try the dav cluster (see below)
set dav=0

if ($dav) then 
    # Couldn't compile hwrf_gettrk properly on dav cluster - Jul 13, 2018
    module load slurm
    cat <<END > $batchfile
#!/bin/tcsh
#SBATCH -A NMMM0013
#SBATCH -t $walltime
#SBATCH -J mpas_gettrk

setenv TMPDIR /glade/scratch/$USER/temp
mkdir -p $TMPDIR

~ahijevyc/bin/mpas_ll_GRIB1.csh $*
END
    sbatch $batchfile
else



    cat <<END > $batchfile
#!/bin/tcsh
#PBS -A NMMM0013
#PBS -N mpas_gettrk
#PBS -q share
#PBS -l walltime=$walltime
#PBS -l select=1:mpiprocs=1:mem=10GB
#PBS -k eod
#PBS -e $TMPDIR/
#PBS -o $TMPDIR/

setenv MPI_USE_ARRAY false

~ahijevyc/bin/mpas_ll_gettrk.csh $*
END
    /opt/pbs/default/bin/qsub $batchfile

endif
