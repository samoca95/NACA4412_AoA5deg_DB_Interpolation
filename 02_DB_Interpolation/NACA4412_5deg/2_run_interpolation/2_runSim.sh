#!/bin/bash

#SBATCH --account=upv103       # Account to launch the job
#SBATCH --job-name=Interp      # Job name
#SBATCH --qos=gp_resa          # Slurm Qos name (associated to partition)
#SBATCH --chdir=.	       # Working directory

#SBATCH --ntasks=112             # Number of tasks
##SBATCH --ntasks-per-node=112  # 112 cpus x node (max. 200 nodes @ gp_resa)
#SBATCH --cpus-per-task=2      # Threads per task (Ncores = ntasks x cpus-per-task)
##SBATCH --hint=nomultithread   # Disallow multithreading to save CPU time
##SBATCH --threads-per-core=1
##SBATCH --exclusive            # Exclusive (automatic if N>1)
##SBATCH --nodes=1              # Number of nodes (max. 200 nodes @ gp_resa)
#SBATCH --time=10:00:00        # Max execution time (DD-HH:MM:SS) (max. 72h @ gp_resa)

#SBATCH --output=nek5000.%j.out
#SBATCH --error=nek5000.%j.err
##SBATCH --mail-type=none       #{begin|end|all|none}
##SBATCH --mail-user=user@mail.com

set -ex

module purge

module load cmake/3.29.2
module load gcc/13.2.0

module load ucx/1.16.0-gcc
module load openmpi/5.0.5-gcc
module load metis/5.1.0-gcc
module load parmetis/4.0.3-gcc-ompi

#source bsc_project upv103

# Path SetUp
#------------------------------------------------
SOURCE_ROOT="$(pwd)"
RUN_PATH="$SOURCE_ROOT"
#------------------------------------------------

cd "$RUN_PATH"
echo "Located at: $RUN_PATH"
echo "$(date +"%d %b %Y at %H:%M") INFO: Start Simulation"
touch SESSION.NAME
echo "naca_wing" > SESSION.NAME
pwd              >> SESSION.NAME

rm *.sch* -f
rm logRun.out -f

# [!IMPORTANT]: If running with srun
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
# [!IMPORTANT]: If running with mpirun
# export SLURM_CPU_BIND=none

srun -n 112 ./nek5000 >> logRun.out 2>&1

echo "$(date +"%d %b %Y at %H:%M") INFO: Simulation End"

# Clean-up the main path
#------------------
# mv *slurm* history/
#-------------------
# bash cpdata.sh

