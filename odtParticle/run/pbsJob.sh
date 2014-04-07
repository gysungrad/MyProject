#!/bin/bash

# This is a simple template pbs submission script for one machine

#PBS -l nodes=16:ppn=8,pmem=300mb,walltime=04:00:00

#PBS -l qos=dol4
#PBS -N myODTjobName

cd $PBS_O_WORKDIR

###################################################

echo "the start time is"
date

mkdir RUNTIME/caseA
mkdir ../input/caseA
#cp ../input/a.inp ../input/odtParam.inp
mkdir ../data/caseA
mpirun -np 128 odt.x
mv ../data/data_* ../data/caseA
mv RUNTIME/runtime_* RUNTIME/caseA

echo "the end time is"
date


exit 0



