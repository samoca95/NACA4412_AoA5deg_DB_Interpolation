# MATLAB for post-processing the binary file 

`main_post_proc_prof.m`: process the `int_fld` data into `.mat` files and export the quantities we are interested in. 

`export_stateOnly`: If one need incorporate more quantities, just modify the [export_drag_decomposition](./matlab_func/export_drag_decomposition.m) script and run 

      export_stateOnly

The files `stat_case_naca4412_200k` in [outputs_data](../outputs_data/) will be modified. 