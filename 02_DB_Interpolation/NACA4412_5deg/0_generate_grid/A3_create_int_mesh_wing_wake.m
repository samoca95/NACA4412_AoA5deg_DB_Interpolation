%% create the interpolating mesh for profiles 
%% From Marco Atzori

clc
clear
close all
addpath matlab_func/
format long

logfile = 'log_wake.txt';
logID = fopen(logfile,'w');

%% compute wing profile

% profile
m=4;
p=4;
xx=12;

%% compute profiles

Rec = 200e3;        % Chord Reynolds number
nu  = 1/Rec;        % Kinematic viscosity
utz = 0;            % Friction velocity

xstart = 1.0;
xend   = 3.0;
Dx_eta = 9;

ystart = -0.1;
yend   = 0.25;
Dy_eta = 5;

zstart = 0;
zend   = 0.2;
Dz_eta = 3; % Adjusted so that final box lies below Dz_eta = 5 (see boxplot)

% Kolmogorov length scale
data = fileread('eta_x.txt'); 
lines = strsplit(data, '\n');     % Split the file content into lines
line1 = strtrim(lines{1});        % Get the first line and trim whitespace
line1 = extractAfter(line1, ':'); % Discard everything up to the separator ":"
eta_x_val = sscanf(line1, '%f;'); % Extract float values separated by ";"
line2 = strtrim(lines{2});        % Get the second line and trim whitespace
line2 = extractAfter(line2, ':'); % Discard everything up to the separator ":"
x_eta = sscanf(line2, '%f;');     % Extract float values separated by ";"

data = fileread('eta_y.txt'); 
lines = strsplit(data, '\n');     % Split the file content into lines
line1 = strtrim(lines{1});        % Get the first line and trim whitespace
line1 = extractAfter(line1, ':'); % Discard everything up to the separator ":"
eta_y_val = sscanf(line1, '%f;'); % Extract float values separated by ";"
line2 = strtrim(lines{2});        % Get the second line and trim whitespace
line2 = extractAfter(line2, ':'); % Discard everything up to the separator ":"
y_eta = sscanf(line2, '%f;');     % Extract float values separated by ";"

figure
hold on
plot(x_eta, eta_x_val, 'k', 'LineWidth', 2)
hold off
print('eta_x.jpg','-dpng','-r300')

figure
hold on
plot(eta_y_val, y_eta, 'k', 'LineWidth', 2)
plot([2.5e-3, 2.5e-3], [y_eta(1), y_eta(end)], 'r--', 'LineWidth', 2)
hold off
print('eta_y.jpg','-dpng','-r300')

% Generate x vector based on eta evolution
xvec(1) = xstart;
xlast = xvec(1);
while xlast < xend
    local_eta_x  = interp1(x_eta, eta_x_val, xlast);
    local_dx     = Dx_eta * local_eta_x;
    xlast        = xvec(end) + local_dx;
    xvec(end+1) = xlast;
end

% Generate y vector based on eta evolution
yvec(1) = ystart;
ylast = yvec(1);
while ylast < yend
    local_eta_y  = interp1(y_eta, eta_y_val, ylast);
    if local_eta_y > 2.5e-3
        local_eta_y = 2.5e-3;
    end
    local_dy     = Dy_eta * local_eta_y;
    ylast        = yvec(end) + local_dy;
    yvec(end+1) = ylast;
end

% Generate z vector
dz    = Dz_eta * mean(eta_x_val);
nz    = ceil((zend-zstart)/dz);
zvec  = linspace(zstart, zend, nz);

figure
hold on
scatter(xvec,ones(length(xvec)), 'k', 'LineWidth', 2)
hold off
print('x_vec.jpg','-dpng','-r300')

figure
hold on
scatter(ones(length(yvec)),yvec, 'k', 'LineWidth', 2)
hold off
print('y_vec.jpg','-dpng','-r300')

nx = length(xvec);
ny = length(yvec);
% nz = length(zvec);

n_points = 0;
x_pts = zeros(nx*ny,1);
y_pts = zeros(nx*ny,1);
% Store the points in the mesh
for ix = 1:nx
    for iy = 1:ny
        n_points = n_points + 1;
        x_pts(n_points,1) = xvec(ix);
        y_pts(n_points,1) = yvec(iy);
    end
end

mesh_wake.Rec = Rec;
mesh_wake.utz = utz;
mesh_wake.x_pts = x_pts;
mesh_wake.y_pts = y_pts;
% Extrude the mesh in the spanwise direction
npts = length(x_pts); % Points in the mesh
mesh_wake.x_pts = repmat(mesh_wake.x_pts, nz, 1);
mesh_wake.y_pts = repmat(mesh_wake.y_pts, nz, 1);
mesh_wake.z_pts = repelem(zvec(:), npts);
mesh_wake.alpha = zeros(size(xvec)); 
mesh_wake.xs = xvec;
mesh_wake.yn = yvec;
mesh_wake.xc = xvec;
mesh_wake.z  = zvec;

Nx = length(mesh_wake.xc); 
Ny = length(mesh_wake.yn);
Nz = length(mesh_wake.z);

disp(['Ntotal=', num2str(length(mesh_wake.x_pts))]);
disp(['Nx=', num2str(Nx), ' Ny=', num2str(Ny), ' Nz=', num2str(Nz)]);
disp(['x/c range: ', num2str(min(mesh_wake.xc)), ' - ', num2str(max(mesh_wake.xc))]);
disp(['y/c range: ', num2str(min(mesh_wake.yn)), ' - ', num2str(max(mesh_wake.yn))]);
disp(['z/c range: ', num2str(min(mesh_wake.z)), ' - ', num2str(max(mesh_wake.z))]);
disp(['xmin=', num2str(min(mesh_wake.x_pts)), ' xmax=', num2str(max(mesh_wake.x_pts))]);
disp(['ymin=', num2str(min(mesh_wake.y_pts)), ' ymax=', num2str(max(mesh_wake.y_pts))]);
disp(['zmin=', num2str(min(mesh_wake.z_pts)), ' zmax=', num2str(max(mesh_wake.z_pts))]);

fprintf(logID,'Ntotal=%d\n',length(mesh_wake.x_pts));
fprintf(logID,'Nx=%d Ny=%d Nz=%d\n',Nx,Ny,Nz);
fprintf(logID,'\n');
fprintf(logID,'x/c range: %f - %f\n', min(mesh_wake.xc), max(mesh_wake.xc));
fprintf(logID,'y/c range: %f - %f\n', min(mesh_wake.yn), max(mesh_wake.yn));
fprintf(logID,'z/c range: %f - %f\n', min(mesh_wake.z), max(mesh_wake.z));
fprintf(logID,'xmin=%f xmax=%f\n', min(mesh_wake.x_pts), max(mesh_wake.x_pts));
fprintf(logID,'ymin=%f ymax=%f\n', min(mesh_wake.y_pts), max(mesh_wake.y_pts));
fprintf(logID,'zmin=%f zmax=%f\n', min(mesh_wake.z_pts), max(mesh_wake.z_pts));
fprintf(logID,'\n');

%% export mesh

int_mesh=mesh_wake;
save('int_mesh_wake.mat',"int_mesh",'-mat')

save_inter_interpolation_mesh('./ZSTRUCT/',mesh_wake);

%%

figure
[xu_plot, yu_plot, ~, xl_plot, yl_plot, ~, ~] = naca_prof(m, p, xx, 100, 1e-6);
% h(1) = scatter(mesh_suction_side.x_pts(plot_mask),mesh_suction_side.y_pts(plot_mask),'.'); hold on
h(1) = scatter(mesh_wake.x_pts,mesh_wake.y_pts,'.'); hold on
plot(xu_plot, yu_plot, 'k','LineWidth',1); hold on
plot(xl_plot, yl_plot, 'k','LineWidth',1); hold on
xlabel('x/c'); ylabel('y/c');
title('Interpolating mesh for the wing profiles')
legend('Wake')
axis equal
hold off
saveas(gcf,'int_mesh_wake.png');

figure
h(1) = scatter(mesh_wake.x_pts,mesh_wake.y_pts,'.'); hold on
xlabel('x/c'); ylabel('y/c');
title('Interpolating mesh for the wing profiles')
legend('Wake')
xlim([1 1.5])
ylim([-0.1 0.25])
axis equal
hold off
saveas(gcf,'int_mesh_wake_detailed.png');


disp('done')
fprintf(logID,'created: %s\n', datetime("now"));
fclose(logID);

%% plot h/eta

hx = xvec(2:end) - xvec(1:end-1);
hy = yvec(2:end) - yvec(1:end-1);
hz = zvec(2) - zvec(1); % Constant in z direction
hx = repmat(hx, ny-1, 1);
hy = repmat(hy, nx-1, 1)';
hz = hz .* ones(ny-1, nx-1);

center_x = (xvec(2:end) + xvec(1:end-1)) / 2;
center_y = (yvec(2:end) + yvec(1:end-1)) / 2;

% Find h/eta for each point
load('eta.mat'); % eta, XC, YN from a thiner interpolation
[XQ, YQ] = meshgrid(center_x, center_y);
local_eta = interp2(XC, YN, real(eta)', XQ, YQ, 'nearest'); % Original grid is much finer

hx_eta = hx ./ local_eta;
hy_eta = hy ./ local_eta;
hz_eta = hz ./ local_eta;
h_eta  = (hx_eta .* hy_eta .* hz_eta) .^ (1/3);

hx_eta(isnan(hx_eta)) = 0;
hy_eta(isnan(hy_eta)) = 0;
hz_eta(isnan(hz_eta)) = 0;
h_eta(isnan(h_eta))   = 0;

figure
hold on
category = repelem({'hx/eta', 'hy/eta', 'hz/eta'}, numel(hx_eta))'; % Repeat each category for all values
values = [reshape(hx_eta, 1, []), reshape(hy_eta, 1, []), reshape(hz_eta, 1, [])]'; % Flatten and concatenate values
tmp = table(categorical(category), values, 'VariableNames', {'Category', 'Values'}); % Convert category to categorical
boxchart(tmp.Category,tmp.Values,'BoxFaceColor',"#0072BD")
yline(9,'--k','h/eta = 9')
yline(5,'--k','h/eta = 5')
hold off
ylabel('h/eta')
set(gca,'FontSize',16)
set(gca,'FontName','Arial')
print('h_boxplot.jpg', '-djpeg', '-r300');
