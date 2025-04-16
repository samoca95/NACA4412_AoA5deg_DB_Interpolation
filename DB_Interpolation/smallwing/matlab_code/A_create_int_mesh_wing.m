%% create the interpolating mesh for profiles

clc
clear
close all
addpath matlab_func/
format long

%% compute wing profile

% profile
m=4;
p=4;
xx=12;

%% compute profiles

Rec=400e3; nu=1/Rec; utz=0.0471; dy1=0.25*nu/utz; %Spacing of first point
                
yn=0.1; %Thicknesses of wall-normal profile delta1

r=1.05; 
% r=1.0025; 
% r=1.002; 

nle=512; 
% nle=1024; 
% nle=2048; 
% nle=4096; 
        
xstart=0.4;
xend=0.8;

tol=1e-6;
[x_pts,y_pts,npoints,yn,alphau,alphal,xc] = create_profiles_wing_interval(m,p,xx,nle,tol,dy1,yn,r,xstart,xend);

%

z=0.01;

mask_suction_side=(x_pts>0)&(y_pts>0);
mask_pressure_side=(x_pts>0)&(y_pts<0);

mesh_suction_side.x_pts=x_pts(mask_suction_side);
mesh_suction_side.y_pts=y_pts(mask_suction_side);
mesh_suction_side.z_pts=ones(size(mesh_suction_side.x_pts))*z;
mesh_suction_side.alpha=alphau;
mesh_suction_side.yn=yn;
mesh_suction_side.xc=xc;
mesh_suction_side.z=z;

mesh_pressure_side.x_pts=x_pts(mask_pressure_side);
mesh_pressure_side.y_pts=y_pts(mask_pressure_side);
mesh_pressure_side.z_pts=ones(size(mesh_pressure_side.x_pts))*z;
mesh_pressure_side.alpha=alphal;
mesh_pressure_side.yn=yn;
mesh_pressure_side.xc=xc;
mesh_pressure_side.z=z;

%%
figure
scatter(mesh_suction_side.x_pts,mesh_suction_side.y_pts,'o'); axis equal; hold on
scatter(mesh_pressure_side.x_pts,mesh_pressure_side.y_pts,'o'); axis equal; hold on

%% export mesh

% save_inter_stat_mesh(path,npoints,x_pts,y_pts);
save_inter_interpolation_mesh('../run_interpolation/ZSTRUCT/',mesh_suction_side);

%%

int_mesh=mesh_suction_side;
% int_mesh=mesh_pressure_side;

save('int_mesh.mat',"int_mesh",'-mat')

disp('done')




