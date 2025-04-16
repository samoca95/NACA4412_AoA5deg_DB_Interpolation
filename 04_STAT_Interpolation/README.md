# Post-processing statistics of NACA4412 at $Re = 200,000$ and $AoA=5^{\circ}$

This repository is basically a fork from [Yuning's repo](https://github.com/Fantasy98/Post-STAT-NACA4412-AoA5deg) with some minor modifications. It contains the post-processing of the simulation results using *Nek5000*. The most of the code are from Ricardo and Marco. 
The framework comprises two ends: 
1. **NEK5000**: **(1) On the fly: in-situ spanwise averaging** and **(2) Post-processing: Ensemble averaging**
2. **Matlab**: **(1) Interpolation map** and **(2) Rotation and Calculation of the raw data** 

## Get Started: Pipeline of the framework: 

1. **Generate the interpolation map use matlab:**

    Generate the mesh with the routines of the `02_DB_Interpolation/NACA4412_5deg/0_generate_grid/` folder, as those ones are more updated and allow for more uniform dx spacing.
    
    Run the MATLAB script `B_create_stat_mesh_wing.m`. Then:

    1. Copy the `{x|y}.fort` files into `01_NEK/ZSTAT`
    2. Copy the `int_stat_mesh.mat` to the `00_Generate_grid` folder for later use.
    

2. Run NEK for **Ensemble average** and **Spectral Interpolation**.

    1. Firstly, collect all the `stsnaca_wing` files and sort them in order. This is done in another folder (`03_STAT_Gathering`)

    2. Run the scrip to create links to the original sts files so that they start by 1:

        bash 00-cp-sts-files.sh

    3. Compile the code and prepare the parameters

            ./01-compile-command naca_wing

      **IMPORTANT** Please Modify the `naca_wing.rea` `PARAM068`, `PARAM069`, `PARAM070` according to your statistics files. 
    
    4. Launch the post-processing code: 

            ./02-run-postNek naca_wing

    Now there will be a binary file `INTERP/int_fld` which contains the quantities and derivatives we need. 

3. **[SIMPLE ALTERNATIVE] Obtain the mean velocity fields and hdf5 files from the `int_fld` file.**

    Inside `02_Convert_hdf5`, run:

        python main_post_proc_prof.py

    This script will output the averaged velocity components and pressure field to an hdf5 file inside `hdf5_files/`.

    Set the parameter `case_type` to "none" to avoid projection of the vector fields to the wall refererence frame (this is done, as in the snapshots interpolation, to givemore flexibility). 
   

4. **[ADVANCED ALTERNATIVE] Use the MATLAB for process the `int_fld` file.**

    (From Yuning, not adapted)

    Launch your MATLAB inside `02_Matlab_code/`, run: 

            main_post_proc_prof 

    The distribution of skin-friction coefficient $c_f$ will be printed afterwards.
  
    Now you can visualize the results in `03_Visualization/`

        python 01-Integral-Quantities.py 
        python 02-WallNormal-Profiles.py 
        python 03-clcd.py 


## Introduction 
+ During the simulation, the span-wise average of the computational field is calculated. The files with prefix `sts` will be dumped and accumulated. 

+ Each file comprises quantities and its derivatives. To check the information, please do 

      head -n 1 stsnaca_wing0.f00001

    Which should give outputs: 

    `#std 8 12 12  1    4692   4692  0.2617500000059E+02    5000  0  1 XS66`,

    where `XS66` indicates there are **66** quantities inside, accompanying with **68** derivatives. 

+ We use a user-defined map (by MATLAB) to acquire the quantities at the points of interest by spectral interpolation in NEK. 
The post-processing mode of NEK will be used for Ensemble average and Spectral interpolation. 

+ The MATLAB functions provide the way to process the quantities into integral quantities and wall-normal profiles for the user to investigate the physics. 

## Notes:

+ If you encounter any issue with compilation. Please unpack the Solver `Nek1093_dong` by: 

        rm -rf Nek1093_dong && tar -xvf Nek1093_dong.tar.gz