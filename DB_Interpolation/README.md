# Nek5000_Interpolation

> [!IMPORTANT] 
> The interpolated fields are in the absolute XY frame of reference to give more freedom to the end user. The x and y components of the vectors should be projected based on the surface angle of the profile if you want to use a wall-tangent and wall-normal reference frame. A routine for that is found in `./3_convert_to_hdf5/B_read_int_data.py`.

Folder structure:

* DEC_Nek5000: Nek5000 source code
* NACA4412_5deg: interpolation folder for the NACA case
* smallwing: small example for a coarser wing

## DB recopilation

Refer to the specific repo to run the recompilation of the raw solution folders. 

This repo will use the `la2naca_wing0.f*` files from the recompilation folder and create symbolic links to them while renaming to start by `.f00001`. It is a requirement of the code to start by 1.

## Steps for interpolation

1. Create the interpolation grid (inside `./0_generate_grid/`):

    Adjust and run the MATLAB `A_create_int_mesh_wing.m` script. This will create files `{x|y|z}.fort3d` inside `./2_run_interpolation/ZSTRUCT/`
  
2. Compile the interpolator (inside `./1_compile_interpolation/`):
   
    Adjust the `SIZE` file's parameters (lp, lhis, IntSize). The rest of parameters (p order, etc, must to be the same as in the original simulation).
   
    Adjust the `.rea` file's parameters:
       - `p011`: the number of time steps should be 0. We will only use it as post-processing.
       - `p068`: to be equal to the number of snapshots to be read and interpolated (content of the database).
    
    **Important:** If you are going to interpolate the db in batches, each batch needs to be renamed so that its first snapshot is 1.
  
    Compile the code by running: 
    
        bash preSim.sh
    
    This will move the necessary files to the run folder too.
    
3. Run the interpolator (inside `./run_iterpolation/`):
   
    Create and order the links to the `la2naca_wing0.f*` snapshots to be interpolated, run:
    
        bash 1_cp_la2_files.sh

    Run the interpolator, use the \# of processors specified in `lp` inside `SIZE`:

        bash 2_runSim.sh

    Note: you might need to adjust the loaded packages depending on the environment.

    Resulting files:
    - `intHistory.f{%5d}`: explains the distribution of points along the processors used
    - `intV{X|Y|Z}.f{%5d}`: interpolated velocity components
    - `intO{X|Y|Z}.f{%5d}`: interpolated vorticity components (O = Omega)
    - `intLA.f{%5d}`: interpolated lamda2 scalar field

4. Convert the results from the interpolator (`./3_convert_to_hdf5/`):
   
        python B_read_int_data.py

    The converted `.h5` files will be stored in the `./3_convert_to_hdf5/hdf5_files/` directory

There are some loading and visualization routines inside the `4_visualization` folder.



