function [] = save2DIntMeshNek(pathmesh,xx,yy)

x_pts=zeros([length(xx)*length(yy),1]);
y_pts=zeros([length(xx)*length(yy),1]);
npoints=0;
for i=1:length(xx)
    for j=1:length(yy)
        npoints=npoints+1;
        x_pts(npoints,1)=xx(i);
        y_pts(npoints,1)=yy(j);
    end
end

%Save x data points in Fortran binary format
fid=fopen([pathmesh,'/ZSTAT/x.fort'],'w','ieee-le.l64');

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
fid=fopen([pathmesh,'/ZSTAT/y.fort'],'w','ieee-le.l64');

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

end