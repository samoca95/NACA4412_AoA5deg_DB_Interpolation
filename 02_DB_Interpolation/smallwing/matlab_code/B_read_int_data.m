%% read interpolated data

clc
clear
close all
addpath matlab_func/

%% import mesh

load('int_mesh.mat')

%%

Nx=length(int_mesh.xc);
Ny=length(int_mesh.yn);
Nz=length(int_mesh.z);

Npoints=Nx*Ny*Nz;

sinalpha=sin(int_mesh.alpha);
cosalpha=cos(int_mesh.alpha);

PATH_INT_FIELDS='../run_interpolation/';

N_FILES=2;

IF_SUCTION_SIDE=0;
IF_PRESSURE_SIDE=0;
IF_VELOCITY=1;
IF_VORTICITY=0;
IF_LA2=0;

%% LOAD INPUT
disp('LOAD INPUT')
for i_file=1:N_FILES
%Load file containing number of wall faces per core
fHistoryname         = strcat(PATH_INT_FIELDS,'intHistory.f',num2str(i_file,'%12.5d'));
% fHistoryname         = strcat('/cfs/klemming/nobackup/a/atzori/MADRID/tDuct_Test/run_360/intHistory_duct0.f',num2str(i_file,'%12.5d'));
history = load(fHistoryname);
ncores=history(end,1)+1; % number of cores

disp(ncores)
disp(i_file)

if IF_VELOCITY==1
variablename='intVX.f';
bb_buff_vx=load_int_file(PATH_INT_FIELDS, variablename, i_file, history, ncores);

variablename='intVY.f';
bb_buff_vy=load_int_file(PATH_INT_FIELDS, variablename, i_file, history, ncores);

variablename='intVZ.f';
bb_buff_vz=load_int_file(PATH_INT_FIELDS, variablename, i_file, history, ncores);
end

if IF_LA2==1
variablename='intLA.f';
bb_buff_la=load_int_file(PATH_INT_FIELDS, variablename, i_file, history, ncores);
end

if IF_VORTICITY==1
variablename='intOX.f';
bb_buff_ox=load_int_file(PATH_INT_FIELDS, variablename, i_file, history, ncores);

variablename='intOY.f';
bb_buff_oy=load_int_file(PATH_INT_FIELDS, variablename, i_file, history, ncores);

variablename='intOZ.f';
bb_buff_oz=load_int_file(PATH_INT_FIELDS, variablename, i_file, history, ncores);
end

%% ROTATE

if IF_VELOCITY==1
test_vx=reshape(bb_buff_vx,Nx,Ny,Nz);
test_vy=reshape(bb_buff_vy,Nx,Ny,Nz);
end

if IF_VORTICITY==1
test_ox=reshape(bb_buff_ox,Nx,Ny,Nz);
test_oy=reshape(bb_buff_oy,Nx,Ny,Nz);
end

if IF_SUCTION_SIDE==1
    for i_k=1:Nz
        for i_y=1:Ny
            for i_x=1:Nx

                if IF_VELOCITY==1
                    vx=test_vx(i_x,i_y,i_k);
                    vy=test_vy(i_x,i_y,i_k);
                    test_vx(i_x,i_y,i_k)=vx*cosalpha(i_x)+vy*sinalpha(i_x);
                    test_vy(i_x,i_y,i_k)=vy*cosalpha(i_x)-vx*sinalpha(i_x);
                end

                if IF_VORTICITY==1
                    ox=test_ox(i_x,i_y,i_k);
                    oy=test_oy(i_x,i_y,i_k);
                    test_ox(i_x,i_y,i_k)=ox*cosalpha(i_x)+oy*sinalpha(i_x);
                    test_oy(i_x,i_y,i_k)=oy*cosalpha(i_x)-ox*sinalpha(i_x);
                end

            end
        end
    end
end

if IF_PRESSURE_SIDE==1
    for i_k=1:Nz
        for i_y=1:Ny
            for i_x=1:Nx


                if IF_VELOCITY==1
                    vx=test_vx(i_x,i_y,i_k);
                    vy=test_vy(i_x,i_y,i_k);
                    test_vx(i_x,i_y,i_k)=vx*cosalpha(i_x)+vy*sinalpha(i_x);
                    test_vy(i_x,i_y,i_k)=vx*sinalpha(i_x)-vy*cosalpha(i_x);
                end

                if IF_VORTICITY==1
                    vx=test_vx(i_x,i_y,i_k);
                    vy=test_vy(i_x,i_y,i_k);
                    test_vx(i_x,i_y,i_k)=vx*cosalpha(i_x)+vy*sinalpha(i_x);
                    test_vy(i_x,i_y,i_k)=vx*sinalpha(i_x)-vy*cosalpha(i_x);
                end

            end
        end
    end
end

if IF_VELOCITY==1
    bb_buff_vx=reshape(test_vx,Nx*Ny*Nz,1,1);
    bb_buff_vy=reshape(test_vy,Nx*Ny*Nz,1,1);
end

if IF_VORTICITY==1
    bb_buff_ox=reshape(test_ox,Nx*Ny*Nz,1,1);
    bb_buff_oy=reshape(test_oy,Nx*Ny*Nz,1,1);
end
end

test_vx=reshape(bb_buff_vx,Ny,Nx,Nz);
test_vy=reshape(bb_buff_vy,Ny,Nx,Nz);
test_vz=reshape(bb_buff_vz,Ny,Nx,Nz);

[XX,YY]=meshgrid(int_mesh.xc,int_mesh.yn);

%%

figure
scatter(int_mesh.x_pts,int_mesh.y_pts); hold on; axis tight equal

%%

set(0,'defaultTextInterpreter','latex'); %trying to set the default
lwidth=4; % line width
legend_on=0;  % show the legend
title_on=0;   % show the title
linewidth_axis=2;  % width of the axis line
scatter_with=1;    % width of the scatter line
figx=1200; % figure horizontal size
figy=500; % figure vertical size
font_size=28;      % font size
set(0, 'DefaultAxesFontWeight', 'normal', ...
    'DefaultAxesFontSize', font_size, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'normal', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'normal', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;

figure('rend','painters','pos',[10 10 figx figy]);
h=pcolor(XX,YY,test_vx(:,:));
set(h, 'EdgeAlpha', 0.1);
axis equal tight
set(gca,'TickLabelInterpreter', 'latex');
title('$u_t=U_t+u_t''$','Interpreter','latex')
xlabel('$x/c$','Interpreter','latex')
ylabel('$y/c$','Interpreter','latex')

figure('rend','painters','pos',[10 10 figx figy]);
h=pcolor(XX,YY,test_vy(:,:));
set(h, 'EdgeAlpha', 0.1);
axis equal tight
set(gca,'TickLabelInterpreter', 'latex');
title('$v_n=V_n+v_n''$','Interpreter','latex')
xlabel('$x/c$','Interpreter','latex')
ylabel('$y/c$','Interpreter','latex')

figure('rend','painters','pos',[10 10 figx figy]);
h=pcolor(XX,YY,test_vz(:,:));
set(h, 'EdgeAlpha', 0.1);
axis equal tight
set(gca,'TickLabelInterpreter', 'latex');
title('$w=W+w''$','Interpreter','latex')
xlabel('$x/c$','Interpreter','latex')
ylabel('$y/c$','Interpreter','latex')

