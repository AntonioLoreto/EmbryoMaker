#!/bin/bash
#SBATCH --job-name="12EMD"
#SBATCH --account=project_2003732
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hugo.cano@uab.cat
#SBATCH --partition=large              
#SBATCH --time=48:00:00          
#SBATCH --error=evo.err        
#SBATCH --output=evo.out    
#SBATCH --nodes=5              
#SBATCH --ntasks=200             # <<< number of parallel tasks --> you must reuse this number for GREASY_NWORKERS  
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G

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

parfile=$(readlink -f $1)
workdir=$(readlink -f $2)
runtime=$3
substitutions=$4

##functions##
zero_pad_u2 () {
    printf "%05d" $1
}

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
   module load hyperqueue

   # Set the directory which hyperqueue will use
   export HQ_SERVER_DIR=$PWD/hq-server-$SLURM_JOB_ID
   mkdir -p "$HQ_SERVER_DIR"

   echo "STARTING HQ SERVER, log in $HQ_SERVER_DIR/HQ.log"
   echo "===================="
   hq server start &>> "$HQ_SERVER_DIR/HQ.log" &
   until hq job list &>/dev/null ; do sleep 1 ; done

   echo "STARTING HQ WORKERS ON $SLURM_NNODES nodes"
   echo "===================="
   srun --exact --cpu-bind=none --mpi=none hq worker start --cpus=${SLURM_CPUS_PER_TASK} &
   hq worker wait "${SLURM_NTASKS}"

   #srun --cpu-bind=none --mpi=none hq worker start --cpus=$SLURM_CPUS_PER_TASK &>> "$HQ_SERVER_DIR/HQ.log" &

   #until [[ $num_up -eq 200 ]]; do
    #   num_up=$(hq worker list 2>/dev/null | grep -c RUNNING )
    #   echo "WAITING FOR WORKERS TO START ( $num_up / 200 )"
    #   sleep 1
   #done

   ## Here you run your submit commands, workflow managers etc...
   ## hq submit <hq submit args> --cpus <n> <COMMAND/executable> <args to program>
   echo "SUBMITING ARRAY JOB"
   pwd
else
    ntasks=${USER_NTASKS:-6}
fi



cp $rang "used_ranges${today}.dat" ##HC 6-10-2021


# This function generates the command line arguments for each individual (running )
# dexe.sh #HC 11-9-2020
generate_args () {  
    cd $running  
    c=1                                                                                    #HC 11-9-2020
    for d in */; do                                                                                  #HC 11-9-2020
            if [ -d "$d" ]; then                                                                     #HC 11-9-2020
                cp ${parfile} ${running}/$(zero_pad_u2 $c)/evaparams.dat
                echo "$runall ${running}/$(zero_pad_u2 $c)/ $BINDIR $runtime $(zero_pad_u2 $c)"      #HC 1-12-2020
                (( c++ ))                                                                            #HC 11-9-2020
            fi                                                                                       #HC 11-9-2020
    done                                                                                             #HC 11-9-2020
    echo "$substitute $BINDIR $substitutions "
}

##HC 11-9-2020 mutate in parallel                             #HC 1-12-2020
echo "Starting simulation"
echo "$(generate_args)" > orders.txt                          #HC 1-12-2020
hq submit --log=logfile.log --each-line=orders.txt run.sh

#for f in /scratch/project_2003732/MorphoDym_3/x*
#  do
#    echo "sublist" $f
#    hq submit --log=logfile.log --cpus=all --wait --array=1-200  $f
#  done

while hq job list --all | grep -q "RUNNING\|PENDING"; do
    echo "WAITING FOR JOBS TO FINISH"
    # Adjust the timing here if you get to much output in the Slurm log file
    # Now set to 30 seconds
    sleep 60
done

echo "===================="
echo "DONE"
echo "===================="
echo "SHUTTING DOWN HYPERQUEUE"
echo "===================="
hq worker stop all
hq server stop

rm finished/*fin* 
echo "Simulation Finished"


