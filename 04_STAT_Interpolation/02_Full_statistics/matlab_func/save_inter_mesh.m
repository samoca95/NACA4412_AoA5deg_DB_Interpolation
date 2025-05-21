function [] = save_inter_mesh(path,npoints,x_pts,y_pts)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Save x data points in Fortran binary format
fid=fopen([path,'ZSTAT/x.fort'],'w','ieee-le.l64');

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

%Save y data points in Fortran binary format
fid=fopen([path,'ZSTAT/y.fort'],'w','ieee-le.l64');

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

%% only for instantaneous measure of Cf, Cd, Cl

% %Save y data points in Fortran binary format
% fid=fopen([path,'ZSTAT/z.fort'],'w','ieee-le.l64');
% 
% %First write 4 bytes integer
% data=npoints;
% eor=length(data)*4;
% count=fwrite(fid,eor,'int32');
% count=fwrite(fid,data,'int32');
% count=fwrite(fid,eor,'int32');
% 
% z_pts=ones(size(y_pts))*0.05;
% %Then write npoints reals
% data=z_pts;
% eor=length(data)*8;
% count=fwrite(fid,eor,'int32');
% count=fwrite(fid,data,'float64');
% count=fwrite(fid,eor,'int32');
% fclose(fid);




end

