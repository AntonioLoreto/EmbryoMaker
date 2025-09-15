#!/bin/bash

# Run single individual
#
# Usage: bash exewrap.sh EXE PARFILE PARENTDIR CHILDDIR

# module purge
# module load ...
# This only mutates individuals
#HC 10-9-2020 This bash was written by Hugo Cano and it:
#HC 10-9-2020 -  calls muta to mutate individuals
#HC 10-9-2020 -  creates twins of each the mutated individual
set -e



childdir=$1                     #executable; elli in my case
bindir=$2                       #HC 1-12-2020 (greasy)
runtime=$3
ind=$4

echo 'runall' $childdir

exe=${bindir}/EMaker
muta=${bindir}/muta.e
fit=${bindir}/fit.e
robust=${bindir}/robust.e
rang=${bindir}/ranges.dat



inviable=${childdir}/inviable.e 
one=${childdir}individual.dat*.dat              #HC 17-11-2020  
filtered=${childdir}/individual.datfitness       #HC 17-11-2020
start=${childdir}/start.dat
OPCval=${childdir}/OPC.val
robval=${childdir}/rob.val
parfile=${childdir}/evaparams.dat

finished=${bindir}/finished
allfin=${finished}/allfin.dat
finish=${finished}/"finish$ind.dat"
START=$(date +%s.%N)

zero_pad_u2 () {
    printf "%05d" $1
}

cd $childdir


while [ ! -f $allfin ]; do
   if [ -f $start ]; then 
      rm $start
      #MUTATION
      timeout 300 $muta individual.dat $parfile $rang || true                            #MUTATION
      if [ ! -f $inviable ]; then
         #DEVELOPMENT
         START=$(date +%s.%N)
         timeout 300000 $exe individual.dat  1 $runtime 1  || true  #HC 10-9-2020        #DEVELOPMENT
         END=$(date +%s.%N)
         DIFF=$(echo "$END - $START" | bc)
         echo $DIFF > "time.dat"   
 
      
         #FITNESS
         if [ ! -f $filtered ]; then                                                      #HC 26-11-2020 Check if the individual was filtered by EMaker or nonmutated
            if [ ! -f $one ]; then                                                        #HC 26-11-2020 Check if there is an output to run fit.e
               echo "OOOPS something went wrong with this individual"                     #HC 26-11-2020
               echo $childir                                                              #HC 30-11-2020
               echo "0.00" > "individual.datfitness"                                      #HC 26-11-2020
               cp ${childdir}individual.dat  ${childdir}individual.err                    #HC 26-11-2020
            else                                                                          #HC 26-11-2020    
               #HC 4-9-2020 Calling the program that calculates robustness                          #HC 9-10-2020
               ## FITNESS
               timeout 3000 $fit $one $parfile $rang || true                              #HC 26-11-2020 Calculate fitness (OPC or EMD)
               ## ROBUSTNESS
               timeout 3000 $robust $one || true #HC 1-9-2020                            #HC 26-11-2020
            fi                                                                            #HC 21-9-2021
         fi                                                                               #HC 26-11-2020

      else
         echo "this individual is inviable"
         echo "0.00" > "individual.datfitness"                                            #HC 26-11-2020
      fi  
                                                                                          
      if [ ! -f $robval ]; then    #HC 26-11-2020
         echo "0.00" > "rob.val"   #HC 26-11-2020
      fi                           #HC 26-11-2020
      
      if [ ! -f $filtered ]; then    #HC 26-11-2020
         echo "0.00" > "individual.datfitness"   #HC 26-11-2020
      fi                           #HC 26-11-2020      


      echo $ind  > $finish   #HC 26-11-2020 this creates the signal file for subs.sh to substitute this individual
      
      START=$(date +%s.%N)
   else
      END=$(date +%s.%N)
      DIFF=$(echo "$END - $START" | bc)
      echo $DIFF > "time_idle.dat"   
   fi
done

echo "You can start" > "start.dat"


