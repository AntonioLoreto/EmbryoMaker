#!/bin/bash

usage="
Usage: bash bootstrap.sh STARTFILE WORKDIR POPSIZE

Sets up the 0th generation,

WORKDIR/
        00000/
                00000/
                    individual.dat
                    individual.datfitness
                00001/
                    individual.dat
                    individual.datfitness
               ...
              POPSIZE-1/
                    individual.dat
                    individual.datfitness
"

set -e

startfile=$(readlink -f $1)
workdir=$(readlink -f $2)
nsize=$3
ncpus=$4
generation=0
seeding=${workdir}/seeding.e
population=${workdir}/population
finished=${workdir}/finished
running=${workdir}/running
zero_pad_u2 () {
    printf "%05d" $1
}

echo $startfile
mkdir ${workdir}/best
mkdir ${workdir}/finished
mkdir ${workdir}/population
mkdir ${workdir}/running

echo "substitution   individual   parent_subs parent   OPC    Symmetry_robustness   mutation   geni   genj   prevalue   newvalue time" > ${workdir}/output_model.dat

for (( i = 1; i <= ncpus; i++ )); do
    individual=${running}/$(zero_pad_u2 $i)
    mkdir $individual
    cp $startfile $individual/individual.dat
    cp $workdir/fourth_original_neighs.dat $individual/fourth_original_neighs.dat   # HC 23-9-2021 this is to calculate joint entropy
    cp $workdir/fourth_original_neighs2.dat $individual/fourth_original_neighs2.dat # HC 23-9-2021
    cp $workdir/ranges.dat $individual/ranges.dat                                   # HC 11-11-2021
    echo "you can start" > $individual/start.dat
    echo "        0         0" > $individual/parent.dat

done

echo $nsize > ${population}/npop.dat
for (( i = 1; i <= nsize; i++ )); do
    c=$(zero_pad_u2 $i)
    cp $startfile ${population}/"individual$c.dat"
    echo "0" >> ${population}/"popids.dat"
   # bc -l <<< "scale=4 ; ${RANDOM}/327670" > ${population}/individual.datfitness
   # cat ${population}/individual.datfitness >> ${population}/population.datfitness
   # rm ${population}/individual.datfitness
done

echo "0" > ${finished}/counter.dat

timeout 3000 $seeding $workdir $startfile $ncpus $nsize

