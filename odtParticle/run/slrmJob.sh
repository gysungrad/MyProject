#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=24:00:00                    # walltime
#SBATCH --ntasks=128                     # number of processor cores (i.e. tasks)
#SBATCH --nodes=16                       # number of nodes
#SBATCH --mem-per-cpu=900M                 # memory per CPU core
#SBATCH -J "SandL"                         # job name
#SBATCH --gid=fslg_crfrg
##SBATCH -p dol4                            #dol4 m7 m6beta
#SBATCH --mail-user=gysungrad@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
#export PBS_QUEUE=batch

###################################################

function rename {
num=$1
for i in `ls RUNTIME/runtime_*`;
do 
    mv $i "RUNTIME/$caseN/runtime_$num"
    num=$[$num+1]
done
num=$1 
for i in `ls -d ../data/data_*`;
do 
    mv $i "../data/$caseN/data_$num"
    num=$[$num+1]
done

}

###################################################
echo "the start time is"
date
###################################################

nRlz=4


#------------ First realization group
#caseN="G10"
#
##cp ../input/Tracer.inp ../input/particle.inp
#mkdir "RUNTIME/$caseN"
#mkdir "../input/$caseN"
#mkdir "../data/$caseN"
#mpiexec -np 128 odt.x
#mv ../data/data_* "../data/$caseN"
#mv RUNTIME/runtime_* "RUNTIME/$caseN"
#
##------------ Other realization groups
#
#shift=128
#it=1
#while [ $it -lt $nRlz ] ; do
#    it=$[$it+1]
#    mpiexec -np 128 odt.x
#    rename $shift
#    shift=$[$shift+128]
#done
#------------ First realization group
#caseN="Re10000_7mm_60umC_1014"
#caseN="SLtracer_1025"
#caseN="test_0220"
#caseN="SL_tracer_betaP0p025_typeC_0221"
#caseN="test_betaP0p07_C6p1Z6p0_typeC_0226"
caseN="homoSL3_040714"

#cp ../input/Hollow.inp ../input/particle.inp
mkdir "RUNTIME/$caseN"
mkdir "../input/$caseN"
mkdir "../data/$caseN"
mpiexec -np 128 odt.x
mv ../data/data_* "../data/$caseN"
mv RUNTIME/runtime_* "RUNTIME/$caseN"

#------------ Other realization groups

shift=128
it=1
while [ $it -lt $nRlz ] ; do
    it=$[$it+1]
    mpiexec -np 128 odt.x
    rename $shift
    shift=$[$shift+128]
done
#------------ First realization group
#caseN="C10_1003"
#
#cp ../input/Corn.inp ../input/particle.inp
#mkdir "RUNTIME/$caseN"
#mkdir "../input/$caseN"
#mkdir "../data/$caseN"
#mpiexec -np 128 odt.x
#mv ../data/data_* "../data/$caseN"
#mv RUNTIME/runtime_* "RUNTIME/$caseN"
#
##------------ Other realization groups
#
#shift=128
#it=1
#while [ $it -lt $nRlz ] ; do
#    it=$[$it+1]
#    mpiexec -np 128 odt.x
#    rename $shift
#    shift=$[$shift+128]
#done
##------------ First realization group
#caseN="S10_1003"
#
#cp ../input/Solid.inp ../input/particle.inp
#mkdir "RUNTIME/$caseN"
#mkdir "../input/$caseN"
#mkdir "../data/$caseN"
#mpiexec -np 128 odt.x
#mv ../data/data_* "../data/$caseN"
#mv RUNTIME/runtime_* "RUNTIME/$caseN"
#
##------------ Other realization groups
#
#shift=128
#it=1
#while [ $it -lt $nRlz ] ; do
#    it=$[$it+1]
#    mpiexec -np 128 odt.x
#    rename $shift
#    shift=$[$shift+128]
#done

#-----------------------------

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
#OUTFILE=""
#mpiexec /fslhome/dol4/compute/particles/C/run/odt.x 

####################################################

echo "the end time is"
date
###################################################

exit 0

# this was the output of the "qsub --show pbsJob.sh" command
#sbatch -N16 -n128 -t01:00:00 --mem-per-cpu=300M --account=dol4 -J "myODTjobName" -o "myODTjobName.o%j" -e "myODTjobName.e%j" 
