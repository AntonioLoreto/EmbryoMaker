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
population=$3
generation=0
seeding=$workdir/seeding.e

zero_pad_u2 () {
    printf "%05d" $1
}

echo $startfile
outdir=${workdir}/$(zero_pad_u2 $generation)
mkdir ${workdir}/best
mkdir -p $outdir
cd $outdir

for (( i = 0; i < population; i++ )); do
    individual=$(zero_pad_u2 $i)
    mkdir $individual
    cp $startfile $individual/individual.dat
    cp $startfile $individual/individual.dat666.dat
    cp $workdir/fourth_original_neighs.dat $individual/fourth_original_neighs.dat   # HC 23-9-2021 this is to calculate joint entropy
    cp $workdir/fourth_original_neighs2.dat $individual/fourth_original_neighs2.dat # HC 23-9-2021
    cp $workdir/ranges.dat $individual/ranges.dat                                   # HC 11-11-2021
    echo "0.00000000000000000" > $individual/OPC.val
    echo "0.00000000000000000  0.00000000000000000" > $individual/rob.val
#    bc -l <<< "scale=4 ; ${RANDOM}/327670" > ${individual}/individual.datfitness
#    echo $i >>${outdir}/ancestors.dat
done

timeout 3000 $seeding $workdir $startfile $population

