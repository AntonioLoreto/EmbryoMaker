#!/bin/bash

# AL 7-3-24>> All comments by Antonio Loreto 

usage="
Usage: bash bootstrap.sh STARTFILE WORKDIR POPSIZE [NCPUS 4 DEVELOPMENT] [RANDOMSEED] [SEEDFILE]

Sets up the 0th generation.

WORKDIR/
        running/
                00000/
                    individual.dat
                    individual.datfitness
                00001/
                    individual.dat
                    individual.datfitness
                            .
                            .
                            .
[NCPUS 4 DEVELOPMENT]/
                    individual.dat
                    individual.datfitness
"

# Exit script immediately if any command fails
set -e

# Validate Arguments
if [ $# -lt 3 ]; then
    echo "$usage"
    exit 1
fi

# Helper function for zero padding
zero_pad_u2() {
    printf "%06d" $1
}

# Read arguments
startfile=$(realpath "$1")  # Ensure correct path resolution
workdir=$(realpath "$2")    # Ensure correct path resolution
nsize=$3                    # Population size
ncpus=${4}              
randomseed=${5:-0}          # Default to 0 if not provided
pathrandseed=${6:-0}        # Default to 0 if not provided

generation=0
seeding="${workdir}/seeding.e"
population="${workdir}/population"
finished="${workdir}/finished"
running="${workdir}/running"

robustatstart="${workdir}/assymetryatstart.e"      # AL 15-4-24>> Robust at start
#complexityatstart="${workdir}/complexityatstart.e"  # AL 15-4-24>> Complexity at start
target="${workdir}/target_1.dat"
developed_ico="${workdir}/chidatttico_evo.dat8209.dat"

# Create necessary directories if they don't exist, and print a warning if they already exist
for dir in "${workdir}/best" "${workdir}/finished" "${workdir}/population" "${workdir}/running"; do
    if [ -d "$dir" ]; then
        echo "[WARNING]: Directory $dir already exists."
    else
        mkdir -p "$dir"
    fi
done

#mkdir "${finished}/all_morfs"           # AL 4-2-25>> Just for testing

echo "    subs    |    individ    |    parent_subs    |    parent    |    fitness    |    bil_sym    |    vol    |      devtime[s]      |      recnum      |      recdonorsubnumber" > "${workdir}/output_model.dat"  # AL 1-4-24

#./fitatstartproj.e icoepi_nuclei_24G.dat10000.dat target_1.dat

./fitatstar.e $developed_ico target_1.dat #AL 19-5-25 
initialfitness=$(cat fitnessatstart.dat)

echo "  0   0   0   0    $initialfitness  0   0   0   0   0   " >> "${workdir}/output_model.dat"

# Set up directories and files for running individuals
for (( i = 1; i <= ncpus; i++ )); do
    individual="${running}/$(zero_pad_u2 $i)"
    mkdir -p "$individual"
    cp "$startfile" "$individual/individual.dat"
    cp "${workdir}/ranges.dat" "$individual/ranges.dat"
    cp "${workdir}/target_1.dat" "$individual/target_1.dat"  # AL 12-11-24: For trait selection
    echo "you can start" > "$individual/start.dat"
    echo "        0         0" > "$individual/parent.dat"
done

# Create population and initialize files
echo "$nsize" > "${population}/npop.dat"
for (( i = 1; i <= nsize; i++ )); do
    c=$(zero_pad_u2 $i)
    cp "$startfile" "${population}/individual$c.dat"
    echo "0" >> "${population}/popids.dat"
done

# AL 15-4-24: Uncomment if necessary
 $robustatstart $developed_ico
# $complexityatstart "$startfile"

#asymmetry and fitness data
for (( i = 1; i <= nsize; i++ )); do
    paste assymetryatstart.dat fitnessatstart.dat -d " " >> fitnassy.txt    # AL 15-4-24
    echo "$initialfitness" >> "${population}/population.datfitness"
done
mv fitnassy.txt "${population}/fitnassy.txt"  # AL 15-4-24: This will be used for inheritance of fitness and asymmetry

# Initialize the counter of substitutions file in the finished directory
echo "0" > "${finished}/counter_subs.dat"

# Initialize the counter of substitutions file in the finished directory
echo "0" > "${finished}/counter_recombinations.dat"

echo -n "" > "${finished}/recombination_history.txt" # AL 26-11-24

# Log path to random seed
echo "path to random seed: $pathrandseed"

# Run the seeding process
timeout 3000 "$seeding" "$workdir" "$startfile" "$ncpus" "$nsize" "$randomseed" "$pathrandseed"
