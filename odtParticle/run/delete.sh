#!/bin/bash

# This is a simple template pbs submission script for one machine

#PBS -l nodes=1:ppn=1,pmem=300mb,walltime=10:00:00

##PBS -l qos=dol4
##PBS -q dol4private
#PBS -N myODTjobName

cd $PBS_O_WORKDIR

###################################################

rm -rf ../data/homoSL* RUNTIME/*
exit 0
