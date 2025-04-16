function [] = save_inter_interpolation_mesh(path,intmesh)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x_pts=intmesh.x_pts;
y_pts=intmesh.y_pts;
z_pts=intmesh.z_pts;

npoints=length(x_pts);

disp(['number of points in the interpolation grid = ',num2str(npoints)])

disp(['xmin=',num2str(min(x_pts)),' xmax=',num2str(max(x_pts))])
disp(['ymin=',num2str(min(y_pts)),' ymax=',num2str(max(y_pts))])
disp(['zmin=',num2str(min(z_pts)),' zmax=',num2str(max(z_pts))])

%% Save x data points in Fortran binary format
fid=fopen([path,'x.fort'],'w','ieee-le.l64');

%First write 4 bytes integer
data=npoints;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=x_pts;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

%% Save y data points in Fortran binary format
fid=fopen([path,'y.fort'],'w','ieee-le.l64');

%First write 4 bytes integer
data=npoints;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=y_pts;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

%% Save z data points in Fortran binary format
fid=fopen([path,'z.fort'],'w','ieee-le.l64');

%First write 4 bytes integer
data=npoints;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=z_pts;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

end

