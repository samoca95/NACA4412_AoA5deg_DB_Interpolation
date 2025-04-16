# Modified
# module swap PrgEnv-cray PrgEnv-intel
echo "Cleaning 01_compile..."
./makenek clean
echo "Compiling case..."
./makenek naca_wing
echo "Done"
