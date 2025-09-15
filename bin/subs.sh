#!/bin/bash


bindir=$1                       #HC 1-12-2020 The bin directory 
substitutions=$2                #HC 1-12-2020 the number of substitutions to be ran
npopulation=$3                  #HC 1-12-2020 the population size substitute #AL 7-3-24>> this argument is actually not passed to the file.

echo 'subs' 
echo $npopulation

substitute=${bindir}/substitute.e #HC 1-12-2020 the exe that makes the substitutions #AL 7-3-24>>this is bad practice because in parallel_evo_cont.sh substitute is = subs.sh which can be confusing.
population=${bindir}/population   #HC 1-12-2020 the directory with the population
finished=${bindir}/finished       #HC 1-12-2020 the directoru with the CPUs that have finished
allfin=${finished}/allfin.dat     #HC 1-12-2020 the file that signals that the simulation is OVER
anyfinished=${finished}/fin*      #HC 1-12-2020 the possible files signaling that a given CPU has finished

cd $finished

sleep 10 #no idea why 10

initial=$(< counter.dat)
final=$(($initial + substitutions))
START=$(date +%s.%N)


while [ ! -f $allfin ]; do                                             #HC 1-12-2020 there are still substitutions to do
       shopt -s nullglob
       files=( fin* )
       if (( ${#files[@]} )); then
         for fin in $finished/fin*; do 
             if [ -f $fin ]; then                                       #HC 1-12-2020 if a CPU has finished #AL 7-3-24>> this means 'if file $fin exists = True'
             sleep 3
             timeout 3000 $substitute $bindir $fin $final $npopulation  #HC 1-12-2020 call the exe that makes the simulations
             START=$(date +%s.%N)
             fi
         done
      else                                                              #HC 1-12-2020 NO CPUs have finished
         sleep 3 #this could be 1      
         END=$(date +%s.%N)
         DIFF=$(echo "$END - $START" | bc)
         echo $DIFF > "time_idle.dat"                                   #HC 1-12-2020 save the time idle
         
      fi
done

