# !/bin/bash

module purge

module load cmake/3.29.2
module load gcc/13.2.0

module load ucx/1.16.0-gcc
module load openmpi/5.0.5-gcc
module load metis/5.1.0-gcc
module load parmetis/4.0.3-gcc-ompi

CASE_NAME='naca_wing'

./compile-script clean
./compile-script ${CASE_NAME}