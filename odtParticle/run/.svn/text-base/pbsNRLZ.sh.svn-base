#!/bin/bash

# This is a simple template pbs submission script for one machine

#PBS -l nodes=16:ppn=8,pmem=300mb,walltime=00:30:00

#PBS -l qos=dol4
##PBS -q dol4private
#PBS -N myODTjobName

cd $PBS_O_WORKDIR

###################################################

caseN="case1"
nRlz=1

echo "the start time is"
date

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

#------------ First realization group

mkdir "RUNTIME/$caseN"
mkdir "../input/$caseN"
mkdir "../data/$caseN"
mpirun -bynode -np 128 odt.x
mv ../data/data_* "../data/$caseN"
mv RUNTIME/runtime_* "RUNTIME/$caseN"

#------------ Other realization groups

shift=128
it=1
while [ $it -lt $nRlz ] ; do
    it=$[$it+1]
    mpirun -np 128 odt.x
    rename $shift
    shift=$[$shift+128]
done

###################################################

echo "the end time is"
date
exit 0



