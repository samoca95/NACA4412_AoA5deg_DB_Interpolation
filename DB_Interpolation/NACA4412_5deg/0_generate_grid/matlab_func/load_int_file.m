function [ bb_buff ] = load_int_file( path_in, variablename, i_file, history, ncores )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fname         = strcat(path_in,variablename,num2str(i_file,'%12.5d'));
[fid,~] = fopen(fname,'r','ieee-le');
hdr           = fread(fid,1,'int32')      ;
F             = fread(fid,hdr,'*char')'   ;
dum5          = fread(fid,1,'*float64')   ;
Rer           = fread(fid,1,'*float64')   ;
Domain        = fread(fid,3,'*float64')   ;
nel           = fread(fid,3,'int32')      ;
Poly          = fread(fid,3,'int32')      ;
% nstat         = fread(fid,1,'int32')      ;
% nderiv        = fread(fid,1,'int32')      ;
time         = fread(fid,1,'*float64')   ;
% timee         = fread(fid,1,'*float64')   ;
% atime         = fread(fid,1,'*float64')   ;
% DT            = fread(fid,1,'*float64')   ;
% nrec          = fread(fid,1,'int32')      ;
% Tint          = fread(fid,1,'*float64')   ;
npoints       = fread(fid,1,'int32')      ;
nfields       = fread(fid,1,'int32')      ;
dum6          = fread(fid,1,'*float64')   ;

disp(fname)
aa=fread(fid,Inf,'*float64'); 
%Assign points taking into account gaps among processors

disp(['read ',fname])
for i=1:ncores
    bb_buff(1+sum(history(1:i-1,2)):sum(history(1:i,2)),1)=aa(1+sum(history(1:i-1,2))+(i-1):sum(history(1:i,2))+(i-1));
    test_read=bb_buff(1+sum(history(1:i-1,2)):sum(history(1:i,2)));
    minrankval=min(test_read);
    maxrankval=max(test_read);
    disp(['irank=',num2str(i-1),' min=',num2str(minrankval),' max=',num2str(max(maxrankval))])
    disp(test_read(1:10))
    disp(test_read(end-10:end))
end


end

