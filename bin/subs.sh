#!/bin/bash
#AL 12-11-2024: This was created by Hugo Cano but modified by Antonio Loreto
bindir=$1                           #HC 1-12-2020 The bin directory 
substitutions=$2                    #HC 1-12-2020 the number of substitutions to be ran

substitute="${bindir}/substitute.e" #HC 1-12-2020 the exe that makes the substitutions
population="${bindir}/population"   #HC 1-12-2020 the directory with the population
finished="${bindir}/finished"       #HC 1-12-2020 the directoru with the CPUs that have finished
allfin="${finished}/allfin.dat"     #HC 1-12-2020 the file that signals that the simulation is OVER
counter_subs="${finished}/counter_subs.dat"
cd "$finished" || { echo "Failed to change directory to $finished"; exit 1; }

echo "subs.sh has been executed. $PWD"
# Read the initial counter
initial=$(< $counter_subs)
final=$((initial + $substitutions))

# Record start time for idle time measurement
START=$(date +%s.%N)

# Loop until allfin is created
while [ ! -f "$allfin" ]; do                                                  #HC 1-12-2020 there are still substitutions to do
    # Find the finished files (sorted by time)
    shopt -s nullglob
    files=(fin*)
    
    if (( ${#files[@]} )); then                                               #HC 1-12-2020 if there are finished files
        for finished_file in $(ls -tr fin*); do                               #AL 5-24-2024 this orders files to be substituted from oldest to newest
            if [ -f "$finished_file" ]; then                                  #HC 1-12-2020 if a CPU has finished
                current=$(< $counter_subs)
                #echo "current substitution: $current"
                #AL 20-11-24:  for recombination we should pass as an argument to substitue.e the probability of recombination
                #that should be read from ranges.dat (this implies modifying the ranges reading) <- 2 be check!!
                timeout 300 "$substitute" "$bindir" "$finished_file" "$final" || { echo "Substitute failed for $finished_file"; continue; }  #HC 1-12-2020 call the exe that makes the simulations
                # Update start time after substitution
                START=$(date +%s.%N)
            fi
        done
    else  #HC 1-12-2020 NO CPUs have finished
        sleep 1
        #echo "sleeping at sub #$current"
        END=$(date +%s.%N)
        DIFF=$(echo "$END - $START" | bc)
        echo "$DIFF" > "time_idle.dat"                                       #HC 1-12-2020 save the time idle 
    fi
done
