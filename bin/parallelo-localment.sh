#!/bin/bash

# Parameter Validation
if [[ -z $1 || -z $2 || -z $3 || -z $4 || -z $5 ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: script.sh [parfile] [BINDIR] [runtime] [substitutions] [totalncpus]"
    exit 1
fi

# Functions 
log() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1"
}

# Parameters and Variables
BINDIR=$(readlink -f $2)
substitute="${BINDIR}/subs.sh"
runall="${BINDIR}/run_all_AL.sh"
running="${BINDIR}/running"
finished="${BINDIR}/finished"
today=$(date +"%Y-%m-%d")   # HC 6-10-2021 Date
rang="${BINDIR}/ranges.dat"
worker="${BINDIR}/worker.sh"
parfile=$(readlink -f "$1") # AL 7-3-24>> first argument is path to paramaters file
runtime=$3
substitutions=$4            # AL 7-3-24>>
totalncpus=$5               # AL 7-3-24>>
population="${BINDIR}/population"

log "Total number of CPUs: $totalncpus"
log "Total number of substitutions: $substitutions"

cp "$rang" "used_ranges${today}.dat" || { echo "Failed to copy ranges file."; exit 1; }

run_task() {
    local index=$1
    local padded_index=$(printf "%06d" $index)
    local run_dir="${running}/${padded_index}/"
    log "Sending to run development CPU #$run_dir"

    mkdir -p "$run_dir" || { echo "Failed to create directory $run_dir"; exit 1; }
    cp "$parfile" "$run_dir/evaparams.dat"
    echo "Let's develop, pal" > "$run_dir/rundevelopment.txt"
    head -n 1 "${population}/fitnassy.txt" > "$run_dir/assynfit.txt"
    bash "$runall" "$run_dir" "$BINDIR" "$runtime" "$padded_index" &
}

# Export necessary variables and function for GNU Parallel
export -f run_task log 
export BINDIR running parfile population runall runtime

# Run tasks in parallel using GNU Parallel
echo "sending development" 
seq 1 $((totalncpus - 1)) | parallel --jobs $((totalncpus - 1)) --joblog joblog.txt run_task {}

echo "sending subs"
nohup bash "$substitute" "$BINDIR" "$substitutions" > output.log 2>&1 &

# Final steps
cd "$BINDIR" || exit
log "Simulation sent."