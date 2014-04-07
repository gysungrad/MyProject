#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=05:00:00                    # walltime
#SBATCH --ntasks=1                     # number of processor cores (i.e. tasks)
#SBATCH --nodes=1                       # number of nodes
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

####################################################
module load matlab
/fslapps/matlab/matlab_7.8/bin/matlab -nodisplay -nojvm -nosplash -r driver 


echo "the end time is"
date
###################################################

exit 0

# this was the output of the "qsub --show pbsJob.sh" command
#sbatch -N16 -n128 -t01:00:00 --mem-per-cpu=300M --account=dol4 -J "myODTjobName" -o "myODTjobName.o%j" -e "myODTjobName.e%j" 
