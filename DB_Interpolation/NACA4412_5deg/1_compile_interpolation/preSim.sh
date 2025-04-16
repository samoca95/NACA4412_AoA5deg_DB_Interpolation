SOURCE_ROOT="$(pwd)"
COMPILE_PATH="$SOURCE_ROOT"
RUN_PATH="$SOURCE_ROOT/../2_run_interpolation"

module purge

module load cmake/3.29.2
module load gcc/13.2.0

module load ucx/1.16.0-gcc
module load openmpi/5.0.5-gcc
module load metis/5.1.0-gcc
module load parmetis/4.0.3-gcc-ompi

cd "$COMPILE_PATH"
bash compile.sh

rm "$RUN_PATH/nek5000" -f
cp nek5000  "$RUN_PATH/nek5000"
echo "INFO: Solver passed!"
cp SIZE naca_wing.rea "$RUN_PATH"
echo "INFO: Setup passed!"
cp naca_wing.wall naca_wing.re2 naca_wing.map  "$RUN_PATH"
echo "INFO: Geometry passed!"

