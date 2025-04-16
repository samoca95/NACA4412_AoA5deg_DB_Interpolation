## Notes

Just follow the numbered scripts 1,2 and 3 to create the symbolic links and check the time separation between snapshots.

Once interpolated the fields, copy them in folder 3_INTERPOLATED_HDF5 and then run script 4 to check that the whole database has the correct DT.

## Transfer from UPV to BSC:

    rsync -avurRtih --progress --stats ./res_1*/la2_data/* oner751450@transfer1.bsc.es:~/scratch/NACA4412_5deg_DB/0_RAW/
