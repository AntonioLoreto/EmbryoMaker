#!/bin/bash
#############################
#SBATCH -J LK2
#SBATCH -n 1
#SBATCH --array=1-160:1
#SBATCH -t 3-00:00:00
#SBATCH -o test_jobsarray-%A-%j-%a.out
#SBATCH -e test_jobsarray-%A-%j-%a.err
#SBATCH -D .
#############################

# This script can be run as a SLURM batch job, or interactively. Default number
# of parallel tasks to run in an interactive job is 2, but can be overridden
# with environment variable USER_NTASKS.

set -e

## Parameters and variables ##
BINDIR=$(readlink -f $2)

substitute=${BINDIR}/subs.sh
runall=${BINDIR}/run_all.sh         ##HC 11-9-2020 This makes mutation   
running=${BINDIR}/running
finished=${BINDIR}/finished
today=$(date +"%Y-%m-%d")     ##HC 6-10-2021 Date
rang=${BINDIR}/ranges.dat
worker=${BINDIR}/worker.sh
parfile=$(readlink -f $1)
workdir=$(readlink -f $2)
runtime=$3
substitutions=$4

##functions##
zero_pad_u2 () {
    printf "%05d" $1
}

fmtID=$(printf "%05d" $SLURM_ARRAY_TASK_ID)
#move files with ancestor trace 
#ancestor_replace () {
#  if (( $1 % 20 == 0 )); then
#    cp $(zero_pad_u2 $1)/ancestors.dat $(zero_pad_u2 $1)/ancestors_gen$1.dat
#    gawk '{print NR-1}' $(zero_pad_u2 $1)/ancestors.dat >$(zero_pad_u2 $1)/temp.dat
#    mv $(zero_pad_u2 $1)/temp.dat $(zero_pad_u2 $1)/ancestors.dat
#  fi
#}

is_slurm_batch_job () {
    [ -n "$SLURM_JOB_ID" ]
}

# start the work. Go to main dir
cd $workdir

# check in what mode the script is running
if is_slurm_batch_job; then
   export SLURM_EXACT=1

   # Load the required modules

   module purge
   module load gnu openmpi/gnu   
   
   echo "LOADED MODULES"
   pwd
else
    ntasks=${USER_NTASKS:-6}
fi



cp $rang "used_ranges${today}.dat" ##HC 6-10-2021


# This function generates the command line arguments for each individual (running )
# dexe.sh #HC 11-9-2020
#generate_args () {  

echo "GENERATING INPUT FILES"
    cd $running  
    c=1                                                                                    #HC 11-9-2020
    for d in */; do                                                                                  #HC 11-9-2020
            if [ -d "$d" ]; then                                                                     #HC 11-9-2020
                cp ${parfile} ${running}/$(zero_pad_u2 $c)/evaparams.dat
                echo " $runall ${running}/$(zero_pad_u2 $c)/ $BINDIR $runtime $(zero_pad_u2 $c)" > $BINDIR/input$(zero_pad_u2 $c).sh     #HC 1-12-2020
                (( c++ ))                                                                            #HC 11-9-2020
            fi                                                                                       #HC 11-9-2020
    done                                                                                             #HC 11-9-2020
    echo "$substitute $BINDIR $substitutions " > $BINDIR/input$(zero_pad_u2 $c).sh
#}

cd $BINDIR
##HC 11-9-2020 mutate in parallel                             #HC 1-12-2020
echo "Starting simulation"
#echo "$(generate_args)" > orders.txt                          #HC 1-12-2020
chmod 777 *sh
srun $worker -i $BINDIR/input$fmtID.sh



echo "===================="
echo "DONE"
echo "===================="

rm finished/*fin* 
echo "Simulation Finished"


