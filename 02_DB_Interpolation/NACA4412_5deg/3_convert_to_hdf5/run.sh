#!/bin/bash

#SBATCH --account=upv103      # Account to launch the job
#SBATCH --job-name=SMC        # Job name
#SBATCH --qos=gp_resa         # Slurm Qos name (associated to partition)
#SBATCH --chdir=.             # Working directory

#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of cores per task
#SBATCH --time=2-24:00:00        # Max execution time (DD-HH:MM:SS) (max. 72h @ gp_resa)

#SBATCH --output=log.%j.out
#SBATCH --error=log.%j.err

set -ex 

script="B_read_int_data.py"

echo "------------------------------------------"
echo "Job started at $(date)"
# # Parse arguments
# while [[ "$#" -gt 0 ]]; do
#     case $1 in
#         --m) script="$2"; shift ;;
#         *) echo "Running script: $1"; exit 1 ;;
#     esac
#     shift
# done
echo "------------------------------------------"

# Load modules
module load anaconda/2024.02

# Run the script
python $script

echo "------------------------------------------"
echo "Job finished at $(date)"
echo "------------------------------------------"