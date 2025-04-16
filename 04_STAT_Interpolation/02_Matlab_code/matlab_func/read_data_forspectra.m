function [top,bottom] = read_data_forspectra(path,x_pts,y_pts,alphal,alphau)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

npoints=length(x_pts);


%Find wall-normal vector and its length
in=find(x_pts==-1);
yn=y_pts(in);
ln=length(yn);

%Number of profiles on top and bottom
np=(npoints/ln)/2;

disp(ln)
disp(np)

%Define structures and coordinates
top=struct;
bottom=struct;


for i=1:np
    top(i).xa=x_pts(ln*(i-1)+1);
    top(i).ya=y_pts(ln*(i-1)+1);
end

bottom(1).xa=x_pts(1);
bottom(1).ya=y_pts(1);

for i=1:np-1
    bottom(i+1).xa=x_pts(ln*(i-1)+1+np*ln);
    bottom(i+1).ya=y_pts(ln*(i-1)+1+np*ln);
end

for i=1:np
    top(i).alpha=alphau(i);
    bottom(i).alpha=alphal(i);
    top(i).yn=yn;
    bottom(i).yn=yn;
end

for i=1:np
    top(i).x=x_pts(ln*(i-1)+1:ln*(i-1)+ln);
    top(i).y=y_pts(ln*(i-1)+1:ln*(i-1)+ln);
end

bottom(1).x=top(1).x;
bottom(1).y=top(1).y;

for i=1:np-1
    bottom(i+1).x=x_pts(ln*(i-1)+1+np*ln:ln*(i-1)+ln+np*ln);
    bottom(i+1).y=y_pts(ln*(i-1)+1+np*ln:ln*(i-1)+ln+np*ln);
end

fname         = [path,'ZSTAT/int_fld'];

[fid,message] = fopen(fname,'r','ieee-le');
hdr           = fread(fid,1,'int32')      ;
F             = fread(fid,hdr,'*char')'   ;
dum5          = fread(fid,1,'*float64')   ;
Rer           = fread(fid,1,'*float64')   ;
Domain        = fread(fid,3,'*float64')   ;
nel           = fread(fid,3,'int32')      ;
Poly          = fread(fid,3,'int32')      ;
nstat         = fread(fid,1,'int32')      ;
nderiv        = fread(fid,1,'int32')      ;
times         = fread(fid,1,'*float64')   ;
timee         = fread(fid,1,'*float64')   ;
atime         = fread(fid,1,'*float64')   ;
DT            = fread(fid,1,'*float64')   ;
nrec          = fread(fid,1,'int32')      ;
Tint          = fread(fid,1,'*float64')   ;
npoints       = fread(fid,1,'int32')      ;
dum6          = fread(fid,1,'*float64')   ;



%Define data arrays
Mat=zeros(ln,2*np);
D=F;

% Arrange stat fields in array form
for ii=1:nstat+nderiv

    if ii == 1
        fseek(fid,0,'cof');
    else
        fseek(fid,8,'cof');
    end

    % Store current field as an array
    MatR = fread(fid,npoints,'*float64');

    % Arrange current field in array form
    for i=1:2*np
        Mat(:,i)=MatR(ln*(i-1)+1:ln*(i-1)+ln);
    end

    if ii<=nstat
        % Generata variable name for field X as FX
        v=genvarname('F', who);
    else
        % Generata variable name for field X as DX
        v=genvarname('D', who);
    end

    % Store matrix in variable FX or Dx
    evalc([v '=INOUT(Mat)']);

end

fclose(fid);



% Simulation parameters
% Bulk Reynolds number Reb=Ub*c/nu
Reb=Rer;

% Domain dimensions
Lx=Domain(1);
Ly=Domain(2);

% Fluid density
rho=1;

% Kinematic viscosity. Both Ub and c are unit normalizing parameters
nu=1/Reb;

% Molecular visosity
mu=rho*nu;



% Fields in the binary record from stat files

%Statistic fields
% u,v,w,p are instantaneous quantities
% averaged in time and homogeneous direction z

% 1.  <u>           % F1
% 2.  <v>           % F2
% 3.  <w>           % F3

% 5.  <uu>          % F5
% 6.  <vv>          % F6
% 7.  <ww>          % F7

% 9.  <uv>          % F9
% 10. <vw>          % F10
% 11. <uw>          % F11

%Derivative fields
% 1. dU/dx           % D1
% 2. dU/dy           % D2
% 3. dV/dx           % D3
% 4. dV/dy           % D4

disp('P')
%Mean pressure. Scalar (MA)
for j=1:2*np-1
    if j<=np
        for i=1:ln
            top(j).mu=mu;
            top(j).nu=nu;
            top(j).rho=rho;            
        end
    else
        for i=1:ln
            bottom(j-np+1).mu=mu;
            bottom(j-np+1).nu=nu;
            bottom(j-np+1).rho=rho;            
        end
        disp(j-np+1)
    end
end

for j=1
    for i=1:ln
        bottom(j).mu=mu;
        bottom(j).nu=nu;
        bottom(j).rho=rho;
    end
end

%Mean velocities. Tensors of Rank 1.
disp('U V W')
for j=1:2*np-1
    if j<=np
        Ralphap=[cos(top(j).alpha) sin(top(j).alpha) 0; ...
            -sin(top(j).alpha) cos(top(j).alpha) 0; ...
            0 0 1];
        for i=1:ln
            prod=Ralphap*[F1(i,j); F2(i,j); F3(i,j)];
            top(j).U(i,1)=prod(1);
            top(j).V(i,1)=prod(2);
            top(j).W(i,1)=prod(3);
        end
    else
         Ralphapp=[cos(bottom(j-np+1).alpha) sin(bottom(j-np+1).alpha) 0; ...
            sin(bottom(j-np+1).alpha) -cos(bottom(j-np+1).alpha) 0; ...
            0 0 1];
        for i=1:ln
            prod=Ralphapp*[F1(i,j); F2(i,j); F3(i,j)];
            bottom(j-np+1).U(i,1)=prod(1);
            bottom(j-np+1).V(i,1)=prod(2);
            bottom(j-np+1).W(i,1)=prod(3);
        end
    end
end

for j=1
    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
        0 0 1];
        for i=1:ln
            prod=Ralphapp*[F1(i,j); F2(i,j); F3(i,j)];
            bottom(j).U(i,1)=prod(1);
            bottom(j).V(i,1)=prod(2);
            bottom(j).W(i,1)=prod(3);
        end
end

disp('u_i u_j')
%Reynolds stress tensor. Tensors of Rank 2.
for j=1:2*np-1
    if j<=np
        Ralphap=[cos(top(j).alpha) sin(top(j).alpha) 0; ...
            -sin(top(j).alpha) cos(top(j).alpha) 0; ...
            0 0 1];
        for i=1:ln
            prod=Ralphap*[(F5(i,j)-F1(i,j).*F1(i,j)) ...
                          (F9(i,j)-F1(i,j).*F2(i,j)) ...
                          (F11(i,j)-F1(i,j).*F3(i,j));
                          (F9(i,j)-F1(i,j).*F2(i,j)) ...
                          (F6(i,j)-F2(i,j).*F2(i,j)) ...
                          (F10(i,j)-F2(i,j).*F3(i,j));
                          (F11(i,j)-F1(i,j).*F3(i,j)) ...
                          (F10(i,j)-F2(i,j).*F3(i,j)) ...
                          (F7(i,j)-F3(i,j).*F3(i,j))]*Ralphap';
            top(j).uu(i,1)=prod(1,1);
            top(j).vv(i,1)=prod(2,2);
            top(j).ww(i,1)=prod(3,3);
            top(j).uv(i,1)=prod(1,2);
            top(j).uw(i,1)=prod(1,3);
            top(j).vw(i,1)=prod(2,3);
        end
    else
         Ralphapp=[cos(bottom(j-np+1).alpha) sin(bottom(j-np+1).alpha) 0; ...
            sin(bottom(j-np+1).alpha) -cos(bottom(j-np+1).alpha) 0; ...
            0 0 1];
        for i=1:ln
            prod=Ralphapp*[(F5(i,j)-F1(i,j).*F1(i,j)) ...
                           (F9(i,j)-F1(i,j).*F2(i,j)) ...
                           (F11(i,j)-F1(i,j).*F3(i,j));
                           (F9(i,j)-F1(i,j).*F2(i,j)) ...
                           (F6(i,j)-F2(i,j).*F2(i,j)) ...
                           (F10(i,j)-F2(i,j).*F3(i,j));
                           (F11(i,j)-F1(i,j).*F3(i,j)) ...
                           (F10(i,j)-F2(i,j).*F3(i,j)) ...
                           (F7(i,j)-F3(i,j).*F3(i,j))]*Ralphapp';
            bottom(j-np+1).uu(i,1)=prod(1,1);
            bottom(j-np+1).vv(i,1)=prod(2,2);
            bottom(j-np+1).ww(i,1)=prod(3,3);
            bottom(j-np+1).uv(i,1)=prod(1,2);
            bottom(j-np+1).uw(i,1)=prod(1,3);
            bottom(j-np+1).vw(i,1)=prod(2,3);
        end
    end
end

for j=1
    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
        0 0 1];
        for i=1:ln
            prod=Ralphapp*[(F5(i,j)-F1(i,j).*F1(i,j)) ...
                           (F9(i,j)-F1(i,j).*F2(i,j)) ...
                           (F11(i,j)-F1(i,j).*F3(i,j));
                           (F9(i,j)-F1(i,j).*F2(i,j)) ...
                           (F6(i,j)-F2(i,j).*F2(i,j)) ...
                           (F10(i,j)-F2(i,j).*F3(i,j));
                           (F11(i,j)-F1(i,j).*F3(i,j)) ...
                           (F10(i,j)-F2(i,j).*F3(i,j)) ...
                           (F7(i,j)-F3(i,j).*F3(i,j))]*Ralphapp';
            bottom(j).uu(i,1)=prod(1,1);
            bottom(j).vv(i,1)=prod(2,2);
            bottom(j).ww(i,1)=prod(3,3);
            bottom(j).uv(i,1)=prod(1,2);
            bottom(j).uw(i,1)=prod(1,3);
            bottom(j).vw(i,1)=prod(2,3);
        end
end



%Velocity gradient tensor. Tensor of Rank 2.
dUdx=D1;
dVdx=D3;
dWdx=D5;

dUdy=D2;
dVdy=D4;
dWdy=D6;

dUdz=zeros(size(dUdx));
dVdz=zeros(size(dUdx));
dWdz=zeros(size(dUdx));

%%

disp('d U_i d x_j')

for j=1:2*np-1
%     Ralphap=eye(3);
%     Ralphapp=eye(3);    
    if j<=np        
        Ralphap=[cos(top(j).alpha) sin(top(j).alpha) 0; ...
            -sin(top(j).alpha) cos(top(j).alpha) 0; ...
            0 0 1];
        for i=1:ln
            prod=Ralphap*[dUdx(i,j) dUdy(i,j) dUdz(i,j);
                          dVdx(i,j) dVdy(i,j) dVdz(i,j);
                          dWdx(i,j) dWdy(i,j) dWdz(i,j)]*Ralphap';

            top(j).dUdx(i,1)=prod(1,1);
            top(j).dUdy(i,1)=prod(1,2);
            top(j).dVdx(i,1)=prod(2,1);
            top(j).dVdy(i,1)=prod(2,2);
            top(j).dWdx(i,1)=prod(3,1);
            top(j).dWdy(i,1)=prod(3,2);
        end
    else
         Ralphapp=[cos(bottom(j-np+1).alpha) sin(bottom(j-np+1).alpha) 0; ...
            sin(bottom(j-np+1).alpha) -cos(bottom(j-np+1).alpha) 0; ...
            0 0 1];
        for i=1:ln
            prod=Ralphapp*[dUdx(i,j) dUdy(i,j) dUdz(i,j);
                           dVdx(i,j) dVdy(i,j) dVdz(i,j);
                           dWdx(i,j) dWdy(i,j) dWdz(i,j)]*Ralphapp';


            bottom(j-np+1).dUdx(i,1)=prod(1,1);
            bottom(j-np+1).dUdy(i,1)=prod(1,2);
            bottom(j-np+1).dVdx(i,1)=prod(2,1);
            bottom(j-np+1).dVdy(i,1)=prod(2,2);
            bottom(j-np+1).dWdx(i,1)=prod(3,1);
            bottom(j-np+1).dWdy(i,1)=prod(3,2);
        end
    end
end

for j=1
    Ralphapp=[cos(bottom(j).alpha) sin(bottom(j).alpha) 0; ...
        sin(bottom(j).alpha) -cos(bottom(j).alpha) 0; ...
        0 0 1];
        for i=1:ln
            prod=Ralphapp*[dUdx(i,j) dUdy(i,j) dUdz(i,j);
                           dVdx(i,j) dVdy(i,j) dVdz(i,j);
                           dWdx(i,j) dWdy(i,j) dWdz(i,j)]*Ralphapp';


            bottom(j).dUdx(i,1)=prod(1,1);
            bottom(j).dUdy(i,1)=prod(1,2);
            bottom(j).dVdx(i,1)=prod(2,1);
            bottom(j).dVdy(i,1)=prod(2,2);
            bottom(j).dWdx(i,1)=prod(3,1);
            bottom(j).dWdy(i,1)=prod(3,2);
        end
end


%%


end

