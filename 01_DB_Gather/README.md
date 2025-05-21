## Notes

Just follow the numbered scripts 1, 2 and 3 to create the symbolic links and check the time separation between snapshots.

Once interpolated the fields, copy them in folder `3_INTERPOLATED_HDF5/` and then run script 4 to check that the whole database has the correct DT.

Note: bash scripts 5 to 7 are created in case you need to join raw interpolation files (if you interpolate from Nek in several batches and you wamt to merge them before converting to python). The folders are joint using symbolic links to avoid memmory issues.