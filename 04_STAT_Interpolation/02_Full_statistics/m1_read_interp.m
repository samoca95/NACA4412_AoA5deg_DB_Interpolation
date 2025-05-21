%% post-processing interpolated data

clc
clear
close all
addpath matlab_func/

%% File names and case

case_name = 'suction'; %'suction'; %'pressure'; %'wake'; %'both';

disp("#######################")
disp("Start Compute Statistics")
disp("#######################")
disp(case_name)

path_binary  = '../01_NEK/INTERP';
binary_file  = 'int_fld';
switch case_name
    case 'suction'
        matfile_name    = './output/naca4412_200k_suction.mat';
    case 'pressure'
        matfile_name    = './output/naca4412_200k_pressure.mat';
    case 'wake'
        matfile_name    = './output/naca4412_200k_wake.mat';
    case 'both'
        matfile_name    = './output/naca4412_200k_both.mat';
end

%% Load the grid used for interpolation
loaded = load('../00_Generate_grid/int_stat_mesh_suction.mat');
mesh = loaded.int_mesh;

%% Read the statistics
long=1;
switch case_name
    case 'both'
        x_pts = mesh.x_pts;
        y_pts = mesh.y_pts;
        xc    = mesh.xc;
        yn    = mesh.yn;
        alphal = mesh.alphal;
        alphal = zeros(size(alphal)); % overwrite alphal with zeros to avoid rotation
        alphau = mesh.alphau;
        alphau = zeros(size(alphau)); % overwrite alphau with zeros to avoid rotation
        [top,bottom] = read_data_profiles(path_binary,binary_file,x_pts,y_pts,alphal,alphau,long,xc);
    otherwise
        x_pts = mesh.x_pts;
        y_pts = mesh.y_pts;
        xc    = mesh.xc;
        yn    = mesh.yn;
        alpha = mesh.alpha;
        alpha = zeros(size(alpha));
        out = read_data_profiles_oneside(case_name,path_binary,binary_file,x_pts,y_pts,alpha,long,xc,yn);
end

%% Save the data
switch case_name
    case 'both'
        save(matfile_name,'top','bottom');
    otherwise
        save(matfile_name,'out');
end
disp('saving done')
disp('done')