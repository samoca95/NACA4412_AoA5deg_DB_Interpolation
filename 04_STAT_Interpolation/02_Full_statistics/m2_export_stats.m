%% post-processing interpolated data

clc
clear
close all
addpath matlab_func/

%% File names and case

case_name = 'suction'; %'suction'; %'pressure'; %'wake'; %'both';

disp("#######################")
disp("Start Exporting Stats")
disp("#######################")
disp(case_name)

path_binary  = '../01_NEK/INTERP';
binary_file  = 'int_fld';
switch case_name
    case 'suction'
        matfile_name    = './output/naca4412_200k_suction.mat';
        output_name     = './output/stat_naca4412_200k_suction.mat';
    case 'pressure'
        matfile_name    = './output/naca4412_200k_pressure.mat';
        output_name     = './output/stat_naca4412_200k_pressure.mat';
    case 'wake'
        matfile_name    = './output/naca4412_200k_wake.mat';
        output_name     = './output/stat_naca4412_200k_wake.mat';
    case 'both'
        matfile_name    = './output/naca4412_200k_both.mat';
        output_name     = './output/stat_naca4412_200k_both.mat';
end

%% Load the grid used for interpolation
loaded = load('../00_Generate_grid/int_stat_mesh_suction.mat');
mesh = loaded.int_mesh;

%% Read the statistics
long = 1;
xc   = mesh.xc;
yn   = mesh.yn;

istart = 1;
iend   = length(xc);
switch case_name
    case 'wake'
        out = load_prof_data_wake(matfile_name, istart, iend, 1);
        disp('Saving stats')
        save(output_name,'out')
    case {'suction', 'pressure'}
        out = load_prof_data_oneside(case_name, matfile_name, istart, iend, 1);
        disp('Saving stats')
        save(output_name,'out')
    case 'both'
        [top, bot] = load_prof_data(matfile_name, istart, iend, 1);
        disp('Saving stats')
        save(output_name,'top','bot')
end
disp('done')

%%
disp("#######################")
disp("     Visualisation     ")
disp("#######################")


% load(output_name)
for i=1:length(xc)
    for j=1:length(yn)
        Dk(i,j) = out(i).Dk(j);
        U(i,j)  = out(i).U(j);
        V(i,j)  = out(i).V(j);
        W(i,j)  = out(i).W(j);
    end
end
figure
hold on
[C,h] = contourf(xc,yn,U',36);
set(h,'LineColor','none')
colorbar()
axis equal
hold off
print('U.jpg', '-djpeg', '-r300');

figure
hold on
[C,h] = contourf(xc,yn,V',36);
set(h,'LineColor','none')
colorbar()
axis equal
hold off
print('V.jpg', '-djpeg', '-r300');

figure
hold on
[C,h] = contourf(xc,yn,W',36);
set(h,'LineColor','none')
colorbar()
axis equal
hold off
print('W.jpg', '-djpeg', '-r300');

figure
hold on
[C,h] = contourf(xc,yn,Dk',36);
set(h,'LineColor','none')
colorbar()
axis equal
hold off
print('Dk.jpg', '-djpeg', '-r300');

disp(['Mean Dk = ', num2str(mean(Dk,[1,2]))])

epsil = -mean(Dk,[1,2]);
nu    = 1/200E3;
eta   = (nu^3/epsil)^0.25;
disp(['eta = ', num2str(eta)])

Lx = 2; Ly = 0.4; Lz = 0.2;
dx = 9*eta; dy = 5*eta; dz = 5*eta;
nx = round(Lx/dx);
ny = round(Ly/dy);
nz = round(Lz/dz);
disp(['nx (9eta) = ', num2str(nx)])
disp(['ny (5eta) = ', num2str(ny)])
disp(['nz (5eta) = ', num2str(nz)])


%% Plotting skin friction coefficient

xaa     = [out(:).xa];
Cfinfs  = [out(:).Cf];
ID      = (xaa <=0.99) & (xaa>=0);

figure;
plot(xaa(ID), Cfinfs(ID), 'k-', 'LineWidth', 3.5);
hold on;

% Add title and labels
title('C_f (Suction Side)');
xlabel('Chordwise Position (x/c)');
ylabel('Skin Friction Coefficient (C_f)');

% Add legend and format it
legend({'Re = 200,000'}, 'Location', 'Best');

% Customize the appearance of the plot
grid on;
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'TickDir', 'out');
set(gca, 'LineWidth', 1);
set(gcf, 'Position', [100, 100, 800, 600]); % Adjust figure size if needed

% Save the figure as an EPS file for inclusion in the paper
print('Cf.jpg', '-djpeg', '-r300');


%% Plotting friction velocity

xaa    = [out(:).xa];
utau   = [out(:).ut];

figure;
plot(xaa(ID), utau(ID), 'k-', 'LineWidth', 3.5);
hold on;

% Add title and labels
title('u_tau (Suction Side)');
xlabel('Chordwise Position (x/c)');
ylabel('Friction Velocity (u_tau)');

% Add legend and format it
legend({'Re = 200,000'}, 'Location', 'Best');

% Customize the appearance of the plot
grid on;
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'TickDir', 'out');
set(gca, 'LineWidth', 1);
set(gcf, 'Position', [100, 100, 800, 600]); % Adjust figure size if needed

% Save the figure as an EPS file for inclusion in the paper
print('u_tau.jpg', '-djpeg', '-r300');

%% Plotting boundary layer thickness

xaa     = [out(:).xa];
delta99 = [out(:).delta99];

figure;
plot(xaa(ID), delta99(ID), 'k-', 'LineWidth', 3.5);
hold on;

% Add title and labels
title('delta99 (Suction Side)');
xlabel('Chordwise Position (x/c)');
ylabel('Boundary Layer Thickness (delta99)');

% Add legend and format it
legend({'Re = 200,000'}, 'Location', 'Best');

% Customize the appearance of the plot
grid on;
box on;
set(gca, 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'TickDir', 'out');
set(gca, 'LineWidth', 1);
set(gcf, 'Position', [100, 100, 800, 600]); % Adjust figure size if needed

% Save the figure as an EPS file for inclusion in the paper
print('delta99.jpg', '-djpeg', '-r300');
