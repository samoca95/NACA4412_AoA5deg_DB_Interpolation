#!/bin/bash

#SBATCH -J wing200k_KEB
#SBATCH -A SNIC2016-34-10

#SBATCH -t 24:00:00
#SBATCH -N 37
#SBATCH --exclusive


# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=yuningw@kth.se


rm *.sch
module load iimpi
mpirun -n 1024  ./nek5000 >> logInit1.txt 2>&1 
