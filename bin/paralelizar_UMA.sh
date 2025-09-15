#!/usr/bin/env bash
# The name to show in queue lists for this job:
#SBATCH -J test_2_EMEV

# Number of desired cores:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Amount of RAM needed for this job:
#SBATCH --mem=8gb

# The time the job will be running, 10 hours:
#SBATCH --time=72:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1

# If you need nodes with special features uncomment the desired constraint line:
##SBATCH --constraint=bigmem
#SBATCH --constraint=cal

# Set output and error files
#SBATCH --error=test_1_EMEV.%A-%a.err
#SBATCH --output=test_1_EMEV.%A-%a.out

# Make an array job, SLURM_ARRAY_TASK_ID will take values from 1 to 10
#SBATCH --array=1-128

# Si está vacío, debido a no usar un array job, entonces se inicializa a 1.
##if test -z "$SLURM_ARRAY_TASK_ID" ; then
##   SLURM_ARRAY_TASK_ID=1
##fi

# Una variable que cambiará entre 0.01 y 0.1, si se lanza con un arrayjob 
##tiempo=`echo $SLURM_ARRAY_TASK_ID|awk '{print  $1/100}'`

# To load some software (you can show the list with 'module avail'):
##module load openmpi_gcc/4.1.1_gcc7.5.0

# the program to execute with its parameters:
#time ./mi_programa argumentos_${tiempo}.jpg > out_${tiempo}.out

# Parameters and Variables
BINDIR=$(readlink -f $2)
substitute="${BINDIR}/subs.sh"
runall="${BINDIR}/run_all_AL.sh"
running="${BINDIR}/running"
finished="${BINDIR}/finished"
today=$(date +"%Y-%m-%d")   # HC 6-10-2021 Date
rang="${BINDIR}/ranges.dat"
parfile=$(readlink -f "$1") # AL 7-3-24>> first argument is path to paramaters file
runtime=$3
substitutions=$4            # AL 7-3-24>>
totalncpus=$5               # AL 7-3-24>>
population="${BINDIR}/population"

run_task() {
    local index=$1
    local padded_index=$(printf "%06d" $index)
    local run_dir="${running}/${padded_index}/"

    # Execute substitute script only for the first index (index == 1)
    if [[ "$index" -eq 1 ]]; then
        echo "Executing substitute script for index 1..."
        time bash "$substitute" "$BINDIR" "$substitutions" > "$run_dir/substitute_output.log" 2>&1
    else 
        # Create directory if necessary
        mkdir -p "$run_dir" || { echo "Failed to create directory $run_dir"; exit 1; }

        # Copy necessary files
        cp "$parfile" "$run_dir/evaparams.dat"
        echo "Let's develop, pal" > "$run_dir/rundevelopment.txt"
        head -n 1 "${population}/fitnassy.txt" > "$run_dir/assynfit.txt"

        # Execute main task
        time bash "$runall" "$run_dir" "$BINDIR" "$runtime" "$padded_index"
    fi
}

# Run the task with the SLURM array index
time run_task ${SLURM_ARRAY_TASK_ID}


