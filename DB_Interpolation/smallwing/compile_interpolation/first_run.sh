#!/bin/bash -l

## Email address and event notification
#SBATCH --mail-user=atzori@mech.kth.se
#SBATCH --mail-type=ALL

#SBATCH -A 2016-34-10

#SBATCH -J Test

#SBATCH -o logfile.%J.out
#SBATCH -e errfile.%J.err 

#SBATCH -n 64
#SBATCH -N 2
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:02:00

## Load any modules you need
#module swap PrgEnv-cray PrgEnv-intel
module load i-compilers/16.0.3 intelmpi/5.1.3

aprun -n ./nek5000
