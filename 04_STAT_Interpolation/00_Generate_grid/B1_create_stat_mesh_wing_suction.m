%% create the interpolating mesh for profiles 
%% From Marco Atzori

clc
clear
close all
addpath matlab_func/
format long

logfile = 'log_stat_suction.txt';
logID = fopen(logfile,'w');

%% compute wing profile

% profile
m=4;
p=4;
xx=12;

%% compute profiles

Rec = 200e3;        % Chord Reynolds number
nu  = 1/Rec;        % Kinematic viscosity
utz = 0.064;        % Friction velocity at x/c=0.4
           

dy1 = 0.5*nu/utz; % Spacing of first point
dy2 = 9*nu/utz;   % Maximum vertical spacing
yn  = 0.12;       % Thicknesses of wall-normal profile delta1
r   = 1.02;       % Expansion ratio from dy1 to dy2, then constant

xstart = 0.2;
xend   = 0.98;
nle    = 575;

tol = 1e-6;
[x_pts,y_pts,npoints,yn,alphau,alphal,xsu,xsl,xc] = create_profiles_wing_interval(m,p,xx,nle,tol,dy1,dy2,yn,r,xstart,xend);

% zvec = 0.1;

mask_suction_side=(x_pts>0)&(y_pts>0);
mask_pressure_side=(x_pts>0)&(y_pts<0);

mesh_suction_side.Rec   = Rec;
mesh_suction_side.utz   = utz;
mesh_suction_side.x_pts = x_pts(mask_suction_side);
mesh_suction_side.y_pts = y_pts(mask_suction_side);
mesh_suction_side.alpha = alphau;
mesh_suction_side.yn    = yn;    % Vector normall to the wall
mesh_suction_side.xs    = xsu;   % Vector tangent to the wall
mesh_suction_side.xc    = xc;    % Vector of chordwise coordinates

% mesh_pressure_side.Rec   = Rec;
% mesh_pressure_side.utz   = utz;
% mesh_pressure_side.x_pts = x_pts(mask_pressure_side);
% mesh_pressure_side.y_pts = y_pts(mask_pressure_side);
% mesh_pressure_side.alpha = alphal;
% mesh_pressure_side.yn    = yn;
% mesh_pressure_side.xs    = xsl;
% mesh_pressure_side.xc    = xc;

%%
% Info about the mesh
Nx = length(mesh_suction_side.xc); 
Ny = length(mesh_suction_side.yn);
[~, ix04]   = min(abs(mesh_suction_side.xc - 0.4));
dxc_04_plus = (mesh_suction_side.xc(ix04+1) - mesh_suction_side.xc(ix04)) * utz / nu;
dxc_plus = diff(mesh_suction_side.xc) * utz / nu;
dxs_plus = diff(mesh_suction_side.xs) * utz / nu;
dxs_04_plus = dxs_plus(ix04);
dyn_plus = diff(mesh_suction_side.yn) * utz / nu;

disp(['Dxc+=', num2str(dxc_04_plus), ' (min=', num2str(min(dxc_plus)), ', max=', num2str(max(dxc_plus)), ')']);
disp(['Dxs+=', num2str(dxs_04_plus), ' (min=', num2str(min(dxs_plus)), ', max=', num2str(max(dxs_plus)), ')']);
disp(['Dyn+=', num2str(dyn_plus(1)), ' to ', num2str(dyn_plus(end))]);
disp(['Nx=', num2str(Nx), ' Ny=', num2str(Ny)]);
disp(['x/c range: ', num2str(min(mesh_suction_side.xc)), ' - ', num2str(max(mesh_suction_side.xc))]);
disp(['y/c range: ', num2str(min(mesh_suction_side.yn)), ' - ', num2str(max(mesh_suction_side.yn))]);
disp(['xmin=', num2str(min(mesh_suction_side.x_pts)), ' xmax=', num2str(max(mesh_suction_side.x_pts))]);
disp(['ymin=', num2str(min(mesh_suction_side.y_pts)), ' ymax=', num2str(max(mesh_suction_side.y_pts))]);

fprintf(logID,'Dxc+=%f (min=%f, max=%f)\n',dxc_04_plus,min(dxc_plus),max(dxc_plus));
fprintf(logID,'Dxs+=%f (min=%f, max=%f)\n',dxs_04_plus,min(dxs_plus),max(dxs_plus));
fprintf(logID,'Dyn+=%f to %f\n',dyn_plus(1),dyn_plus(end));
fprintf(logID,'\n');
fprintf(logID,'Ntotal=%d\n',length(mesh_suction_side.x_pts));
fprintf(logID,'Nx=%d Ny=%d\n',Nx,Ny);
fprintf(logID,'\n');
fprintf(logID,'x/c range: %f - %f\n', min(mesh_suction_side.xc), max(mesh_suction_side.xc));
fprintf(logID,'y/c range: %f - %f\n', min(mesh_suction_side.yn), max(mesh_suction_side.yn));
fprintf(logID,'xmin=%f xmax=%f\n', min(mesh_suction_side.x_pts), max(mesh_suction_side.x_pts));
fprintf(logID,'ymin=%f ymax=%f\n', min(mesh_suction_side.y_pts), max(mesh_suction_side.y_pts));
fprintf(logID,'\n');

% xmin_plot = 0.4;
% xmax_plot = 0.42;
% ymin_plot = 0;
% ymax_plot = 0.25;
% plot_mask = (mesh_suction_side.x_pts > xmin_plot) & (mesh_suction_side.x_pts < xmax_plot) & ...
%             (mesh_suction_side.y_pts > ymin_plot) & (mesh_suction_side.y_pts < ymax_plot);

figure
[xu_plot, yu_plot, ~, xl_plot, yl_plot, ~, ~] = naca_prof(m, p, xx, 100, tol);
% h(1) = scatter(mesh_suction_side.x_pts(plot_mask),mesh_suction_side.y_pts(plot_mask),'.'); hold on
h(1) = scatter(mesh_suction_side.x_pts,mesh_suction_side.y_pts,'.'); hold on
plot(xu_plot, yu_plot, 'k','LineWidth',1); hold on
plot(xl_plot, yl_plot, 'k','LineWidth',1); hold on
% h(2) = scatter3(mesh_pressure_side.x_pts,mesh_pressure_side.y_pts,'o'); axis equal; hold on
xlabel('x/c'); ylabel('y/c');
title('Interpolating mesh for the wing profiles')
legend('Suction side')
% legend('Suction side','Pressure side')
axis equal
hold off
saveas(gcf,'stat_mesh_suction.png');

%% export mesh

save_inter_stat_mesh('./ZSTAT/',mesh_suction_side);

%%
int_mesh=mesh_suction_side;
% int_mesh=mesh_pressure_side;
save('int_stat_mesh_suction.mat',"int_mesh",'-mat')


disp('done')
fprintf(logID,'created: %s\n', datetime("now"));
fclose(logID);




