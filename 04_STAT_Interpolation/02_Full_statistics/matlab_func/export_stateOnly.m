%% post-processing interpolated data

clc
clear
close all

addpath matlab_func/

%%

path='../01_NEK/'; 

% case_name = "oppo";
case_name = "benchmark";
% case_name = 'blowing';
% case_name = 'suction';

% case_name = 'experiments'

disp(case_name)
%%
long=1;
interpolation_mesh_type='ProfilesPostProcessing'; % and FIK % and Budget

istart =18 ;
iend   =80 ; 

disp("#######################")
disp("Start Exporting Drags")
disp("#######################")

switch case_name
    case 'benchmark'
        matfile_name='../outputs_data/benchmark_naca4412_200k.mat';
        output_name='../outputs_data/stat_benchmark_naca4412_200k.mat';

    case 'blowing'
        matfile_name='../outputs_data/blw_ss_01_naca4412_200k.mat';
        output_name='../outputs_data/stat_blw_ss_01_naca4412_200k.mat';
    
    case 'suction'
        matfile_name='../outputs_data/sct_ss_01_naca4412_200k.mat';
        output_name='../outputs_data/stat_sct_ss_01_naca4412_200k.mat';
    
    case 'bodyforce'
        matfile_name='../outputs_data/body_naca4412_200k.mat';
        output_name='../outputs_data/stat_body_naca4412_200k.mat';

    case 'oppo'
        matfile_name='../outputs_data/oppo_naca4412_200k.mat';
        output_name='../outputs_data/stat_oppo_naca4412_200k.mat';
    case 'experiments'
        matfile_name='../outputs_data/exp_naca4412_200k.mat';
        output_name ='../outputs_data/stat_exp_naca4412_200k.mat';
end 

[top200kref,bot200kref]=load_prof_data(matfile_name, istart, iend, 1 );
% Extract the suction side and pressure side 
Re200k.top.ref=export_drag_decomposition(top200kref, istart, iend);
Re200k.bot.ref=export_drag_decomposition(bot200kref, istart, iend);
save(output_name,'Re200k')


disp(case_name)
disp('done')