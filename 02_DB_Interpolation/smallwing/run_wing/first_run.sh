#!/bin/bash

#SBATCH -J w200k_BESKOW
#SBATCH -A 2018-2-4

#SBATCH -t 16:00:00
#SBATCH -N 128
#SBATCH --exclusive
#SBATCH -n 4096
#SBATCH --ntasks-per-node=32

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=marco.atzori.92@gmail.com


rm *.sch
module swap PrgEnv-cray PrgEnv-intel

FOLDER=Run13Quater

srun -n 4096  ./nek5000 >> logRun14Quater_BES_pts.txt 2>&1 
mkdir $FOLDER
cp stsnaca_wing0.f0000* $FOLDER
cp rs8naca_wing0.*	$FOLDER
cp la2naca_wing0.*	$FOLDER
cp ptsnaca_wing0.f000*	$FOLDER
 


