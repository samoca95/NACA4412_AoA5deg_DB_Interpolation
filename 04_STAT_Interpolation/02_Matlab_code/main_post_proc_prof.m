%% post-processing interpolated data

clc
clear
close all

addpath matlab_func/

%%

path='../01_NEK/'; 

% case_name = "oppo";
% case_name = "benchmark";
case_name = 'benchmark';
% case_name = 'suction';
% case_name = 'bodyforce';
% case_name = 'experiments'

disp("#######################")
disp("Start Compute Statistics")
disp("#######################")
disp(case_name)

switch case_name
    case 'benchmark'
        binary_file='int_fld'
        matfile_name='../outputs_data/benchmark_naca4412_200k.mat';
    case 'blowing'
        binary_file='int_fld'
        matfile_name='../outputs_data/blw_ss_01_naca4412_200k.mat';
    case 'suction'
        matfile_name='../outputs_data/sct_ss_01_naca4412_200k.mat';
    
    case 'bodyforce'
        binary_file='int_fld'
        matfile_name='../outputs_data/body_naca4412_200k.mat';
    case 'oppo'
        binary_file= 'int_fld'
        matfile_name='../outputs_data/oppo_naca4412_200k.mat';
    case 'experiments'
        matfile_name='../outputs_data/exp_naca4412_200k.mat';

end 

%%
long=1;
interpolation_mesh_type='ProfilesPostProcessing'; % and FIK % and Budget

load int_mesh.mat;

%%

switch interpolation_mesh_type
    case 'ProfilesPostProcessing'
        
        [top,bottom] = read_data_profiles(path,binary_file,x_pts,y_pts,alphal,alphau,long,xc);
        
        save(matfile_name,'top','bottom');
        disp('saving done')

    otherwise
        disp('NO MESH TYPE')
        return
end

%%

disp('done')


close all 
% clc 
istart =18 ;
iend   =80 ; 

% 


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

% 
[top200kref,bot200kref]=load_prof_data(matfile_name, istart, iend, 1 );

% Extract the suction side and pressure side 
Re200k.top.ref=export_drag_decomposition(top200kref, istart, iend);
Re200k.bot.ref=export_drag_decomposition(bot200kref, istart, iend);
save(output_name,'Re200k')
disp('done')

%%
disp("#######################")
disp("      Visualisation          ")
disp("#######################")


load(output_name)

xaa  = [Re200k.top.ref(:).xa];
Cfinfs = [Re200k.top.ref(:).Cf];
ID =  (xaa <=0.95) & (xaa>=0);

figure;
plot(xaa, Cfinfs, 'k-', 'LineWidth', 3.5);
hold on;

% Add title and labels
title('C_f (Suction Side)');
xlabel('Chordwise Position (x/c)');
ylabel('Skin Friction Coefficient (C_f)');

% Add legend and format it
% legend({'Re = 200,000'}, 'Location', 'Best');

% Customize the appearance of the plot
grid on;
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'TickDir', 'out');
set(gca, 'LineWidth', 1);
set(gcf, 'Position', [100, 100, 800, 600]); % Adjust figure size if needed

% Save the figure as an EPS file for inclusion in the paper
print('suction_side_cf.jpg', '-djpeg', '-r300');

xaa  = [Re200k.bot.ref(:).xa];
Cfinfs = [Re200k.bot.ref(:).Cfinf];
ID =  (xaa <=0.95) & (xaa>=0);

figure;
plot(xaa, Cfinfs, 'k-', 'LineWidth', 3.5);
hold on;

% Add title and labels
title('C_f (Pressure Side)');
xlabel('Chordwise Position (x/c)');
ylabel('Skin Friction Coefficient (C_f)');

% Add legend and format it
% legend({'Re = 200,000'}, 'Location', 'Best');

% Customize the appearance of the plot
grid on;
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'TickDir', 'out');
set(gca, 'LineWidth', 1);
set(gcf, 'Position', [100, 100, 800, 600]); % Adjust figure size if needed

% Save the figure as an EPS file for inclusion in the paper
print('pressure_side_cf.jpg', '-djpeg', '-r300');