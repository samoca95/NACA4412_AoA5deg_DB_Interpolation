function out = read_data_profiles_oneside(side,path,binary_file,x_pts,y_pts,alpha,long,xc,yn)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%Find wall-normal vector and its length
ln=length(yn);

%Number of profiles on top and bottom
np=length(xc);

disp(['ln : ',num2str(ln)])
disp(['np : ',num2str(np)])

%Define structures and coordinates
out=struct;

for i=1:np
    out(i).xa=x_pts(ln*(i-1)+1);
    out(i).ya=y_pts(ln*(i-1)+1);
    out(i).xc=xc(i);
end

for i=1:np
    out(i).alpha=alpha(i);
    out(i).yn=yn;
end

for i=1:np
    out(i).x=x_pts(ln*(i-1)+1:ln*(i-1)+ln);
    out(i).y=y_pts(ln*(i-1)+1:ln*(i-1)+ln);
end


fname         = [path,'/',binary_file];

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

% Check number of points
if npoints ~= ln*np
    disp('Number of points in the file does not match the number of points in the mesh')
    disp(['Number of points in the file: ',num2str(npoints)])
    disp(['Number of points in the mesh: ',num2str(ln*np)])
end

%Define data arrays
Mat=zeros(ln,np);
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
    for i=1:np
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
% 4.  <p>           % F4

% 5.  <uu>          % F5
% 6.  <vv>          % F6
% 7.  <ww>          % F7
% 8.  <pp>          % F8

% 9.  <uv>          % F9
% 10. <vw>          % F10
% 11. <uw>          % F11

% 12. <pu>          % F12
% 13. <pv>          % F13
% 14. <pw>          % F14

% 15. <pdudx>       % F15
% 16. <pdudy>       % F16
% 17. <pdudz>       % F17

% 18. <pdvdx>       % F18
% 19. <pdvdy>       % F19
% 20. <pdvdz>       % F20

% 21. <pdwdx>       % F21
% 22. <pdwdy>       % F22
% 23. <pdwdz>       % F23

% 24. <uuu>         % F24
% 25. <vvv>         % F25
% 26. <www>         % F26
% 27. <ppp>         % F27

% 28. <uuv>         % F28
% 29. <uuw>         % F29
% 30. <vvu>         % F30
% 31. <vvw>  	    % F31
% 32. <wwu>         % F32
% 33. <wwv>         % F33
% 34. <uvw>         % F34

% 35. <uuuu>        % F35
% 36. <vvvv>        % F36
% 37. <wwww>        % F37
% 38. <pppp>        % F38

% 39. <uuuv>        % F39
% 40. <uuvv>        % F40
% 41. <uvvv> 	    % F41

% 42. e11: <((du/dx)^2+(du/dy)^2+(du/dz)^2)> % F42
% 43. e22: <((dv/dx)^2+(dv/dy)^2+(dv/dz)^2)> % F43
% 44. e33: <((dw/dx)^2+(dw/dy)^2+(dw/dz)^2)> % F44
% 45. e12: <(du/dx*dv/dx+du/dy*dv/dy+du/dz*dv/dz)>  % F45
% 46. e13: <(du/dx*dw/dx+du/dy*dw/dy+du/dz*dw/dz)>  % F46
% 47. e23: <(dv/dx*dw/dx+dv/dy*dw/dy+dv/dz*dw/dz)>  % F47

% 48. <dw/dx*dw/dx> % F48
% 49. <dw/dy*dw/dy> % F49
% 50. <dw/dx*dw/dy> % F50

% 51. <du/dx*du/dx> % F51
% 52. <du/dy*du/dy> % F52
% 53. <du/dx*du/dy> % F53

% 54. <dv/dx*dv/dx> % F54
% 55. <dv/dy*dv/dy> % F55
% 56. <dv/dx*dv/dy> % F56

% 57. <du/dx*dv/dx> % F57
% 68. <du/dy*dv/dy> % F58
% 59. <du/dx*dv/dy> % F59
% 60. <du/dy*dv/dx> % F60

% MA: from the forcing terms:
% 61. <fx> % F61
% 62. <fy> % F62
% 63. <fz> % F63
% 64. <fx*ux> % F64
% 65. <fy*uy> % F65
% 66. <fz*uz> % F66


%Derivative fields
% 1. dU/dx           % D1
% 2. dU/dy           % D2
% 3. dV/dx           % D3
% 4. dV/dy           % D4

% 5. dW/dx           % D5
% 6. dW/dy           % D6
% 7. dP/dx           % D7
% 8. dP/dy           % D8

% 9.  d<uu>/dx       % D9
% 10. d<uu>/dy       % D10
% 11. d<vv>/dx       % D11
% 12. d<vv>/dy       % D12

% 13. d<ww>/dx       % D13
% 14. d<ww>/dy       % D14
% 15. d<pp>/dx       % D15
% 16. d<pp>/dy       % D16

% 17. d<uv>/dx       % D17
% 18. d<uv>/dy       % D18
% 19. d<vw>/dx       % D19
% 20. d<vw>/dy       % D20

% 21. d<uw>/dx       % D21
% 22. d<uw>/dy       % D22
% 23. d<uuu>/dx      % D23
% 24. d<uuu>/dy      % D24

% 25. d<vvv>/dx      % D25
% 26. d<vvv>/dy      % D26
% 27. d<www>/dx      % D27
% 28. d<www>/dy      % D28

% 29. d<ppp>/dx      % D29
% 30. d<ppp>/dy      % D30
% 31. d<uuv>/dx      % D31
% 32. d<uuv>/dy      % D32

% 33. d<uuw>/dx      % D33
% 34. d<uuw>/dy      % D34
% 35. d<vvu>/dx      % D35
% 36. d<vvu>/dy      % D36

% 37. d<vvw>/dx      % D37
% 38. d<vvw>/dy      % D38
% 39. d<wwu>/dx      % D39
% 40. d<wwu>/dy      % D40

% 41. d<wwv>/dx      % D41
% 42. d<wwv>/dy      % D42
% 43. d<uvw>/dx      % D43
% 44. d<uvw>/dy      % D44

% 45. d2U/dx2        % D45
% 46. d2U/dy2        % D46
% 47. d2V/dx2        % D47
% 48. d2V/dy2        % D48

% 49. d2W/dx2        % D49
% 50. d2W/dy2        % D50
% 51. d2<uu>/dx2     % D51
% 52. d2<uu>/dy2     % D52

% 53. d2<vv>/dx2     % D53
% 54. d2<vv>/dy2     % D54
% 55. d2<ww>/dx2     % D55
% 56. d2<ww>/dy2     % D56

% 57. d2<uv>/dx2     % D57
% 58. d2<uv>/dy2     % D58
% 59. d2<uw>/dx2     % D59
% 60. d2<uw>/dy2     % D60

% 61. d2<vw>/dx2     % D61
% 62. d2<vw>/dy2     % D62

% 63. d<pu>/dx       % D63
% 64. d<pu>/dy       % D64
% 65. d<pv>/dx       % D65
% 66. d<pv>/dy       % D66

% 67. d<pw>/dx       % D67
% 68. d<pw>/dy       % D68


% Matrices for rotation
switch side
    case {'suction', 'top', 'upper'}
        for j=1:np
            Ralphap(j,1:3,1:3)=[ cos(out(j).alpha) sin(out(j).alpha) 0; ...
                                -sin(out(j).alpha) cos(out(j).alpha) 0; ...
                                 0                 0                 1];
        end
    case {'pressure', 'bot', 'bottom'}
        for j=1:np
            Ralphap(j,1:3,1:3)=[cos(out(j).alpha)  sin(out(j).alpha) 0; ...
                                sin(out(j).alpha) -cos(out(j).alpha) 0; ...
                                0                  0                 1];
        end
    otherwise
        for j=1:np
            Ralphap(j,1:3,1:3)=[1 0 0; ...
                                0 1 0; ...
                                0 0 1];
        end
end

disp('P')
%Mean pressure. Scalar (MA)
for j=1:np
    for i=1:ln
        out(j).P(i,1)=F4(i,j);
        out(j).mu=mu;
        out(j).nu=nu;
        out(j).rho=rho;            
    end
end

%Mean velocities. Tensors of Rank 1.
disp('U V W')
for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[F1(i,j); F2(i,j); F3(i,j)];
        out(j).U(i,1)=prod(1);
        out(j).V(i,1)=prod(2);
        out(j).W(i,1)=prod(3);
    end
end

%Pressure gradient.  Tensors of Rank 1.
disp('dPdx dPdy dPdz')
%Mean velocities. Tensors of Rank 1.
for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[D7(i,j); D8(i,j); D8(i,j)*0.0];
        out(j).dPdx(i,1)=prod(1);
        out(j).dPdy(i,1)=prod(2);
        out(j).dPdz(i,1)=prod(3);
    end
end


disp('u_i u_j')
%Reynolds stress tensor. Tensors of Rank 2.
for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[(F5(i,j)-F1(i,j).*F1(i,j)) ...
                    (F9(i,j)-F1(i,j).*F2(i,j)) ...
                    (F11(i,j)-F1(i,j).*F3(i,j));
                    (F9(i,j)-F1(i,j).*F2(i,j)) ...
                    (F6(i,j)-F2(i,j).*F2(i,j)) ...
                    (F10(i,j)-F2(i,j).*F3(i,j));
                    (F11(i,j)-F1(i,j).*F3(i,j)) ...
                    (F10(i,j)-F2(i,j).*F3(i,j)) ...
                    (F7(i,j)-F3(i,j).*F3(i,j))]*Ralph';
        out(j).uu(i,1)=prod(1,1);
        out(j).vv(i,1)=prod(2,2);
        out(j).ww(i,1)=prod(3,3);
        out(j).uv(i,1)=prod(1,2);
        out(j).uw(i,1)=prod(1,3);
        out(j).vw(i,1)=prod(2,3);
    end
end

%Mean, RMS, skewness and flatness of pressure
P=F4;
pinf = max(max(P))-0.5;
pp=F8-P.*P;
ppp=F27-3*P.*pp-P.*P.*P;
pppp=F38-4*P.*ppp-6*P.*P.*pp-P.*P.*P.*P;

%Normalize pressure
prms=sqrt(pp);
pskew=ppp./(pp).^(3/2);
pflat=pppp./(pp).^(2);

%%

%Assign RMS of pressure
for j=1:np
    for i=1:ln
        out(j).pinf=pinf;
        out(j).P(i,1)=P(i,j);
        out(j).pp(i,1)=pp(i,j);
        out(j).prms(i,1)=prms(i,j);
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

for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[dUdx(i,j) dUdy(i,j) dUdz(i,j);
                    dVdx(i,j) dVdy(i,j) dVdz(i,j);
                    dWdx(i,j) dWdy(i,j) dWdz(i,j)]*Ralph';

        out(j).dUdx(i,1)=prod(1,1);
        out(j).dUdy(i,1)=prod(1,2);
        out(j).dVdx(i,1)=prod(2,1);
        out(j).dVdy(i,1)=prod(2,2);
        out(j).dWdx(i,1)=prod(3,1);
        out(j).dWdy(i,1)=prod(3,2);
        % check stuff
        out(j).dUdxnotR(i,1)=dUdx(i,j);
        out(j).dVdynotR(i,1)=dVdy(i,j);
    end
end


U=F1;
V=F2;
W=F3;

duudx=D9-2*U.*dUdx;
dvvdx=D11-2*V.*dVdx;
dwwdx=D13-2*W.*dWdx;
duvdx=D17-U.*dVdx-V.*dUdx;
duwdx=D21-U.*dWdx-W.*dUdx;
dvwdx=D19-V.*dWdx-W.*dVdx;

duudy=D10-2*U.*dUdy;
dvvdy=D12-2*V.*dVdy;
dwwdy=D14-2*W.*dWdy;
duvdy=D18-U.*dVdy-V.*dUdy;
duwdy=D22-U.*dWdy-W.*dUdy;
dvwdy=D20-V.*dWdy-W.*dVdy;

duudz=zeros(size(duudx));
dvvdz=zeros(size(duudx));
dwwdz=zeros(size(duudx));
duvdz=zeros(size(duudx));
duwdz=zeros(size(duudx));
dvwdz=zeros(size(duudx));

%%

if long==1

%Production tensor. Tensor of Rank 2.
%Definition of the production tensor assuming fully-developed flow, i.e.,
%d()/dz=0, as a function of x and y.
%Pij=-(<uiuk>*dUj/dxk+<ujuk>*dUi/dxk)
%P11=-2*(<uu>*dU/dx+<uv>*dU/dy)
%P12=-(<uu>*dV/dx +<uv>*dV/dy+<uv>*dU/dx+<vv>*dU/dy)
%P13=-(<uu>*dW/dx +<uv>*dW/dy+<uw>*dU/dx+<vw>*dU/dy)
%P22=-2*(<uv>*dV/dx+<vv>*dV/dy)
%P23=-(<uv>*dW/dx+<vv>*dW/dy+<uw>*dV/dx+<vw>*dV/dy)
%P33=-2*(<uw>*dW/dx+<vw>*dW/dy)


%Reynolds stress tensor. Tensor of Rank 2.
R2_tensor(1:3,1:3,1:ln,1:np)=0;
for j=1:np
    for i=1:ln
        R2_tensor(:,:,i,j)=[(F5(i,j)-F1(i,j).*F1(i,j))  (F9(i,j)-F1(i,j).*F2(i,j))  (F11(i,j)-F1(i,j).*F3(i,j));
                            (F9(i,j)-F1(i,j).*F2(i,j))  (F6(i,j)-F2(i,j).*F2(i,j))  (F10(i,j)-F2(i,j).*F3(i,j));
                            (F11(i,j)-F1(i,j).*F3(i,j)) (F10(i,j)-F2(i,j).*F3(i,j)) (F7(i,j)-F3(i,j).*F3(i,j))];
    end
end

uu=squeeze(squeeze(R2_tensor(1,1,:,:)));
vv=squeeze(squeeze(R2_tensor(2,2,:,:)));
ww=squeeze(squeeze(R2_tensor(3,3,:,:)));
uv=squeeze(squeeze(R2_tensor(1,2,:,:)));
uw=squeeze(squeeze(R2_tensor(1,3,:,:)));
vw=squeeze(squeeze(R2_tensor(2,3,:,:)));
%

for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[-2*(uu(i,j).*dUdx(i,j)+uv(i,j).*dUdy(i,j)) ...
                    -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
                    -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j));
                    -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
                    -2*(uv(i,j).*dVdx(i,j)+vv(i,j).*dVdy(i,j)) ...
                    -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j));
                    -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j)) ...
                    -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j)) ...
                    -2*(uw(i,j).*dWdx(i,j)+vw(i,j).*dWdy(i,j))]*Ralph';

        out(j).Pxx(i,1)=prod(1,1);
        out(j).Pyy(i,1)=prod(2,2);
        out(j).Pzz(i,1)=prod(3,3);
        out(j).Pxy(i,1)=prod(1,2);
        out(j).Pxz(i,1)=prod(1,3);
        out(j).Pyz(i,1)=prod(2,3);
    end
end

%Dissipation tensor. Tensor of Rank 2.
%Definition of the dissipation tensor assuming fully-developed flow, i.e.,
%d()/dz=0, as a function of x and y.
%Dij=-2*nu*<dui/dxk*duj/dxk>

%e11_tot=<((du/dx)^2+(du/dy)^2+(du/dz)^2)> % F42
%e22_tot=<((dv/dx)^2+(dv/dy)^2+(dv/dz)^2)> % F43
%e33_tot=<((dw/dx)^2+(dw/dy)^2+(dw/dz)^2)> % F44
%e12_tot=<(du/dx*dv/dx+du/dy*dv/dy+du/dz*dv/dz)>  % F45
%e13_tot=<(du/dx*dw/dx+du/dy*dw/dy+du/dz*dw/dz)>  % F46
%e23_tot=<(dv/dx*dw/dx+dv/dy*dw/dy+dv/dz*dw/dz)>  % F47

%e11=e11_tot-(dU/dx)^2-(dU/dy)^2
%e22=e22_tot-(dV/dx)^2-(dV/dy)^2
%e33=e33_tot-(dW/dx)^2-(dW/dy)^2
%e12=e12_tot-(dU/dx*dV/x)-(dU/dy*dV/dy)
%e13=e13_tot-(dU/dx*dW/x)-(dU/dy*dW/dy)
%e23=e23_tot-(dV/dx*dW/x)-(dV/dy*dW/dy)

e11(1:ln,1:np)=0;
e22(1:ln,1:np)=0;
e33(1:ln,1:np)=0;
e12(1:ln,1:np)=0;
e13(1:ln,1:np)=0;
e23(1:ln,1:np)=0;
%%

for j=1:np
    for i=1:ln
        dUidxj(:,:,i,j)=[dUdx(i,j) dUdy(i,j) dUdz(i,j);
                         dVdx(i,j) dVdy(i,j) dVdz(i,j);
                         dWdx(i,j) dWdy(i,j) dWdz(i,j)];
    end
end

for j=1:np
    e11(:,j)=F42(:,j)- ...
        (squeeze(dUidxj(1,1,:,j))).^2-(squeeze(dUidxj(1,2,:,j))).^2;
    e22(:,j)=F43(:,j)- ...
        (squeeze(dUidxj(2,1,:,j))).^2-(squeeze(dUidxj(2,2,:,j))).^2;
    e33(:,j)=F44(:,j)- ...
        (squeeze(dUidxj(3,1,:,j))).^2-(squeeze(dUidxj(3,2,:,j))).^2;
    e12(:,j)=F45(:,j)- ...
        squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(2,1,:,j))- ...
        squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(2,2,:,j));
    e13(:,j)=F46(:,j)- ...
        squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(3,1,:,j))- ...
        squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(3,2,:,j));
    e23(:,j)=F47(:,j)- ...
        squeeze(dUidxj(2,1,:,j)).*squeeze(dUidxj(3,1,:,j))- ...
        squeeze(dUidxj(2,2,:,j)).*squeeze(dUidxj(3,2,:,j));
end

%
for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[e11(i,j) e12(i,j) e13(i,j);
                    e12(i,j) e22(i,j) e23(i,j);
                    e13(i,j) e23(i,j) e33(i,j)]*Ralph';

        out(j).Dxx(i,1)=-2*nu*prod(1,1);
        out(j).Dyy(i,1)=-2*nu*prod(2,2);
        out(j).Dzz(i,1)=-2*nu*prod(3,3);
        out(j).Dxy(i,1)=-2*nu*prod(1,2);
        out(j).Dxz(i,1)=-2*nu*prod(1,3);
        out(j).Dyz(i,1)=-2*nu*prod(2,3);

    end
end
%%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean forcing terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean forcing terms:
for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[F61(i,j); F62(i,j); F63(i,j)];
        out(j).Forcingx(i,1)=prod(1);
        out(j).Forcingy(i,1)=prod(2);
        out(j).Forcingz(i,1)=prod(3);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Additional dissipation from the RT filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Dissipation from RT filter. Tensors of Rank 1.
for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[F64(i,j); F65(i,j); F66(i,j)];
        out(j).DRTx(i,1)=prod(1);
        out(j).DRTy(i,1)=prod(2);
        out(j).DRTz(i,1)=prod(3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mean convection tensor. Tensor of Rank 2.
%Definition of the mean convection tensor assuming
%fully-developed flow, i.e., d()/dz=0, as a function of x and y.
%Cij=Uk*d<uiuj>/dxk
%Note that under this definition: Production + Dissipation - Convection ...
%C11=U*d(uu)/dx+V*d(uu)/dy
%C22=U*d(vv)/dx+V*d(vv)/dy
%C33=U*d(ww)/dx+V*d(ww)/dy
%C12=U*d(uv)/dx+V*d(uv)/dy
%C13=U*d(uw)/dx+V*d(uw)/dy
%C23=U*d(vw)/dx+V*d(vw)/dy

for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[U(i,j).*duudx(i,j)+V(i,j).*duudy(i,j) ...
                    U(i,j).*duvdx(i,j)+V(i,j).*duvdy(i,j) ...
                    U(i,j).*duwdx(i,j)+V(i,j).*duwdy(i,j);
                    U(i,j).*duvdx(i,j)+V(i,j).*duvdy(i,j) ...
                    U(i,j).*dvvdx(i,j)+V(i,j).*dvvdy(i,j) ...
                    U(i,j).*dvwdx(i,j)+V(i,j).*dvwdy(i,j);
                    U(i,j).*duwdx(i,j)+V(i,j).*duwdy(i,j) ...
                    U(i,j).*dvwdx(i,j)+V(i,j).*dvwdy(i,j) ...
                    U(i,j).*dwwdx(i,j)+V(i,j).*dwwdy(i,j)]*Ralph';

        out(j).Cxx(i,1)=prod(1,1);
        out(j).Cyy(i,1)=prod(2,2);
        out(j).Czz(i,1)=prod(3,3);
        out(j).Cxy(i,1)=prod(1,2);
        out(j).Cxz(i,1)=prod(1,3);
        out(j).Cyz(i,1)=prod(2,3);

    end
end

%Turbulent transport tensor. Tensor of Rank 2.
%Definition of the turbulent transport tensor assuming fully-developed
%flow, i.e., d()/dz=0, as a function of x and y.
%Tij=-d<uiujuk>/dxk
%T11=-(d<uuuk>)/dxk, K=1,3=-[d<uuu>/dx+d<uuv>/dy+d<uuw>/dz]
%T22=-(d<vvuk>)/dxk, K=1,3=-[d<vvu>/dx+d<vvv>/dy+d<vvw>/dz]
%T33=-(d<wwuk>)/dxk, K=1,3=-[d<wwu>/dx+d<wwv>/dy+d<www>/dz]
%T12=-(d<uvuk>)/dxk, K=1,3=-[d<uvu>/dx+d<uvv>/dy+d<uvw>/dz]
%T13=-(d<uwuk>)/dxk, K=1,3=-[d<uwu>/dx+d<uwv>/dy+d<uww>/dz]
%T23=-(d<vwuk>)/dxk, K=1,3=-[d<vwu>/dx+d<vwv>/dy+d<vww>/dz]

duuudx=D23-3*U.*U.*dUdx-3*(U.*duudx+uu.*dUdx);
dvvudx=D35-2*(V.*duvdx+uv.*dVdx)-(U.*dvvdx+vv.*dUdx) ...
    -(V.*V.*dUdx+2*U.*V.*dVdx);
dwwudx=D39-2*(W.*duwdx+uw.*dWdx)-(U.*dwwdx+ww.*dUdx) ...
    -(W.*W.*dUdx+2*U.*W.*dWdx);
duvudx=D31-2*(U.*duvdx+uv.*dUdx)-(V.*duudx+uu.*dVdx) ...
    -(U.*U.*dVdx+2*U.*V.*dUdx);
duwudx=D33-2*(U.*duwdx+uw.*dUdx)-(W.*duudx+uu.*dWdx) ...
    -(U.*U.*dWdx+2*U.*W.*dUdx);
dvwudx=D43-(U.*dvwdx+vw.*dUdx)-(V.*duwdx+uw.*dVdx)-(W.*duvdx+uv.*dWdx) ...
    -(U.*V.*dWdx+U.*W.*dVdx+V.*W.*dUdx);

duuvdy=D32-2*(U.*duvdy+uv.*dUdy)-(V.*duudy+uu.*dVdy) ...
    -(U.*U.*dVdy+2*U.*V.*dUdy);
dvvvdy=D26-3*(V.*dvvdy+vv.*dVdy)-3*V.*V.*dVdy;
dwwvdy=D42-2*(W.*dvwdy+vw.*dWdy)-(V.*dwwdy+ww.*dVdy) ...
    -(W.*W.*dVdy+2*V.*W.*dWdy);
duvvdy=D36-2*(V.*duvdy+uv.*dVdy)-(U.*dvvdy+vv.*dUdy) ...
    -(V.*V.*dUdy+2*U.*V.*dVdy);
duwvdy=D44-(U.*dvwdy+vw.*dUdy)-(V.*duwdy+uw.*dVdy)-(W.*duvdy+uv.*dWdy) ...
    -(U.*V.*dWdy+U.*W.*dVdy+V.*W.*dUdy);
dvwvdy=D38-2*(V.*dvwdy+vw.*dVdy)-(W.*dvvdy+vv.*dWdy) ...
    -(V.*V.*dWdy+2*V.*W.*dVdy);

duuwdz=zeros(size(duuudx));
dvvwdz=zeros(size(duuudx));
dwwwdz=zeros(size(duuudx));
duvwdz=zeros(size(duuudx));
duwwdz=zeros(size(duuudx));
dvwwdz=zeros(size(duuudx));

%%

for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[-(duuudx(i,j)+duuvdy(i,j)+duuwdz(i,j)) ...
                    -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
                    -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j));
                    -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
                    -(dvvudx(i,j)+dvvvdy(i,j)+dvvwdz(i,j)) ...
                    -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j));
                    -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j)) ...
                    -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j)) ...
                    -(dwwudx(i,j)+dwwvdy(i,j)+dwwwdz(i,j))]*Ralph';

        out(j).Txx(i,1)=prod(1,1);
        out(j).Tyy(i,1)=prod(2,2);
        out(j).Tzz(i,1)=prod(3,3);
        out(j).Txy(i,1)=prod(1,2);
        out(j).Txz(i,1)=prod(1,3);
        out(j).Tyz(i,1)=prod(2,3);

    end
end

%Viscous diffusion tensor. Tensor of Rank 2.
%Definition of the viscous diffusion tensor assuming fully-developed
%flow, i.e., d()/dz=0, as a function of x and y.
%VDij=nu*d2(uiuj)/dxk2
%VD11=nu*(d2(uu)/dx2+d2(uu)/dy2+d2(uu)/dz2)
%VD22=nu*(d2(vv)/dx2+d2(vv)/dy2+d2(vv)/dz2)
%VD33=nu*(d2(ww)/dx2+d2(ww)/dy2+d2(ww)/dz2)
%VD12=nu*(d2(uv)/dx2+d2(uv)/dy2+d2(uv)/dz2)
%VD13=nu*(d2(uw)/dx2+d2(uw)/dy2+d2(uw)/dz2)
%VD23=nu*(d2(vw)/dx2+d2(vw)/dy2+d2(vw)/dz2)

d2uudx2=D51-2*(U.*D45+dUdx.*dUdx);
d2vvdx2=D53-2*(V.*D47+dVdx.*dVdx);
d2wwdx2=D55-2*(W.*D49+dWdx.*dWdx);
d2uvdx2=D57-(V.*D45+U.*D47+2*dUdx.*dVdx);
d2uwdx2=D59-(U.*D49+W.*D45+2*dUdx.*dWdx);
d2vwdx2=D61-(V.*D49+W.*D47+2*dVdx.*dWdx);

d2uudy2=D52-2*(U.*D46+dUdy.*dUdy);
d2vvdy2=D54-2*(V.*D48+dVdy.*dVdy);
d2wwdy2=D56-2*(W.*D50+dWdy.*dWdy);
d2uvdy2=D58-(V.*D46+U.*D48+2*dUdy.*dVdy);
d2uwdy2=D60-(U.*D50+W.*D46+2*dUdy.*dWdy);
d2vwdy2=D62-(V.*D50+W.*D48+2*dVdy.*dWdy);

d2uudz2=zeros(size(d2uudx2));
d2vvdz2=zeros(size(d2uudx2));
d2wwdz2=zeros(size(d2uudx2));
d2uvdz2=zeros(size(d2uudx2));
d2uwdz2=zeros(size(d2uudx2));
d2vwdz2=zeros(size(d2uudx2));

for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[nu*(d2uudx2(i,j)+d2uudy2(i,j)+d2uudz2(i,j)) ...
                    nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
                    nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j));
                    nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
                    nu*(d2vvdx2(i,j)+d2vvdy2(i,j)+d2vvdz2(i,j)) ...
                    nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j));
                    nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j)) ...
                    nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j)) ...
                    nu*(d2wwdx2(i,j)+d2wwdy2(i,j)+d2wwdz2(i,j))]*Ralph';

        out(j).VDxx(i,1)=prod(1,1);
        out(j).VDyy(i,1)=prod(2,2);
        out(j).VDzz(i,1)=prod(3,3);
        out(j).VDxy(i,1)=prod(1,2);
        out(j).VDxz(i,1)=prod(1,3);
        out(j).VDyz(i,1)=prod(2,3);

    end
end

%Velocity-pressure-gradient tensor. Tensor of Rank 2.
%Definition of the velocity-pressure-gradient tensor assuming
%fully-developed flow, i.e., d()/dz=0, as a function of x and y.
%Piij=-1/rho*(<ui*dp/dxj>+<uj*dp/dxi>)
%Pi11=-2/rho*<u*dp/dx>
%Pi22=-2/rho*<v*dp/dy>
%Pi33=-2/rho*<w*dp/dz>
%Pi12=-1/rho*(<u*dp/dy>+<v*dp/dx>)
%Pi13=-1/rho*(<u*dp/dz>+<w*dp/dx>)
%Pi23=-1/rho*(<v*dp/dz>+<w*dp/dy>)

%Now, since we don't compute <ui*dp/dxj>, we use the chain rule to express
%these terms as a function of velocity gradients.
%<ui*dp/dxj>=d(<p*ui>)/dxj-<p*dui/dxj>

%We define the pressure transport and pressure strain tensors.

%Pressure transport tensor. Tensor of Rank 2.
%Definition of the pressure transport tensor assuming
%fully-developed flow, i.e., d()/dz=0, as a function of x and y.
%PTij=-1/rho*(<d(p*ui)/dxj>+<d(p*uj)/dxi>)
%PT11=-2/rho*<d(p*u)/dx>
%PT22=-2/rho*<d(p*v)/dy>
%PT33=-2/rho*<d(p*w)/dz>
%PT12=-1/rho*(<d(p*u)/dy>+<d(p*v)/dx>)
%PT13=-1/rho*(<d(p*u)/dz>+<d(p*w)/dx>)
%PT23=-1/rho*(<d(p*v)/dz>+<d(p*w)/dy>)

dpudx=D63-P.*dUdx-U.*D7;
dpvdx=D65-P.*dVdx-V.*D7;
dpwdx=D67-P.*dWdx-W.*D7;

dpudy=D64-P.*dUdy-U.*D8;
dpvdy=D66-P.*dVdy-V.*D8;
dpwdy=D68-P.*dWdy-W.*D8;

dpudz=zeros(size(dpudx));
dpvdz=zeros(size(dpudx));
dpwdz=zeros(size(dpudx));

for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[(dpudx(i,j)) ...
                    (dpudy(i,j)+dpvdx(i,j)) ...
                    (dpudz(i,j)+dpwdx(i,j));
                    (dpudy(i,j)+dpvdx(i,j)) ...
                    (dpvdy(i,j)) ...
                    (dpvdz(i,j)+dpwdy(i,j));
                    (dpudz(i,j)+dpwdx(i,j)) ...
                    (dpvdz(i,j)+dpwdy(i,j)) ...
                    (dpwdz(i,j))]*Ralph';

        out(j).PTxx(i,1)=prod(1,1);
        out(j).PTyy(i,1)=prod(2,2);
        out(j).PTzz(i,1)=prod(3,3);
        out(j).PTxy(i,1)=prod(1,2);
        out(j).PTxz(i,1)=prod(1,3);
        out(j).PTyz(i,1)=prod(2,3);

    end
end

%Pressure strain tensor. Tensor of Rank 2.
%Definition of the pressure strain tensor assuming
%fully-developed flow, i.e., d()/dz=0, as a function of x and y.
%PSij=-1/rho*(<p*dui/dxj>+<p*duj/dxi>)
%PS11=-2/rho*<p*du/dx>
%PS22=-2/rho*<p*dv/dy>
%PS33=-2/rho*<p*dw/dz>
%PS12=-1/rho*(<p*du/dy>+<p*dv/dx>)
%PS13=-1/rho*(<p*du/dz>+<p*dw/dx>)
%PS23=-1/rho*(<p*dv/dz>+<p*dw/dy>)

%<pdudx>       % F15
%<pdudy>       % F16
%<pdudz>       % F17

%<pdvdx>       % F18
%<pdvdy>       % F19
%<pdvdz>       % F20

%<pdwdx>       % F21
%<pdwdy>       % F22
%<pdwdz>       % F23

pdudx=F15-P.*dUdx;
pdudy=F16-P.*dUdy;
pdudz=F17-P.*dUdz;

pdvdx=F18-P.*dVdx;
pdvdy=F19-P.*dVdy;
pdvdz=F20-P.*dVdz;

pdwdx=F21-P.*dWdx;
pdwdy=F22-P.*dWdy;
pdwdz=F23-P.*dWdz;

for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[(pdudx(i,j)) ...
                    (pdudy(i,j)+pdvdx(i,j)) ...
                    (pdudz(i,j)+pdwdx(i,j));
                    (pdudy(i,j)+pdvdx(i,j)) ...
                    (pdvdy(i,j)) ...
                    (pdvdz(i,j)+pdwdy(i,j));
                    (pdudz(i,j)+pdwdx(i,j)) ...
                    (pdvdz(i,j)+pdwdy(i,j)) ...
                    (pdwdz(i,j))]*Ralph';

        out(j).PSxx(i,1)=prod(1,1);
        out(j).PSyy(i,1)=prod(2,2);
        out(j).PSzz(i,1)=prod(3,3);
        out(j).PSxy(i,1)=prod(1,2);
        out(j).PSxz(i,1)=prod(1,3);
        out(j).PSyz(i,1)=prod(2,3);

    end
end

%%

%Construct velocity-pressure-gradient-tensor
for i=1:length(out)
    out(i).Pixx=-2/rho*(out(i).PTxx-out(i).PSxx);
    out(i).Piyy=-2/rho*(out(i).PTyy-out(i).PSyy);
    out(i).Pizz=-2/rho*(out(i).PTzz-out(i).PSzz);
    out(i).Pixy=-1/rho*(out(i).PTxy-out(i).PSxy);
    out(i).Pixz=-1/rho*(out(i).PTxz-out(i).PSxz);
    out(i).Piyz=-1/rho*(out(i).PTyz-out(i).PSyz);
end

%Budget for each component of the Reynolds stress tensor
%Without mean convection
for i=1:length(out)
    out(i).Sxx=out(i).Pxx+out(i).Dxx+out(i).Txx+out(i).VDxx+out(i).Pixx;
    out(i).Syy=out(i).Pyy+out(i).Dyy+out(i).Tyy+out(i).VDyy+out(i).Piyy;
    out(i).Szz=out(i).Pzz+out(i).Dzz+out(i).Tzz+out(i).VDzz+out(i).Pizz;
    out(i).Sxy=out(i).Pxy+out(i).Dxy+out(i).Txy+out(i).VDxy+out(i).Pixy;
    out(i).Sxz=out(i).Pxz+out(i).Dxz+out(i).Txz+out(i).VDxz+out(i).Pixz;
    out(i).Syz=out(i).Pyz+out(i).Dyz+out(i).Tyz+out(i).VDyz+out(i).Piyz;
end

%With mean convection
for i=1:length(out)
    out(i).Scxx=out(i).Pxx+out(i).Dxx+out(i).Txx+out(i).VDxx+out(i).Pixx-out(i).Cxx;
    out(i).Scyy=out(i).Pyy+out(i).Dyy+out(i).Tyy+out(i).VDyy+out(i).Piyy-out(i).Cyy;
    out(i).Sczz=out(i).Pzz+out(i).Dzz+out(i).Tzz+out(i).VDzz+out(i).Pizz-out(i).Czz;
    out(i).Scxy=out(i).Pxy+out(i).Dxy+out(i).Txy+out(i).VDxy+out(i).Pixy-out(i).Cxy;
    out(i).Scxz=out(i).Pxz+out(i).Dxz+out(i).Txz+out(i).VDxz+out(i).Pixz-out(i).Cxz;
    out(i).Scyz=out(i).Pyz+out(i).Dyz+out(i).Tyz+out(i).VDyz+out(i).Piyz-out(i).Cyz;
end

%Skewness tensor. Tensor of Rank 3.
R3_tensor_tot(1:3,1:3,1:3,1:ln,1:np)=0;
R3_tensor(1:3,1:3,1:3,1:ln,1:np)=0;

%Form of the tensor.
% [ uuu uuv uuw ] [ uuv uvv uvw ] [ uuw uvw uww ]
% [ uuv uvv uvw ] [ uvv vvv vvw ] [ uvw vvw vww ]
% [ uuw uvw uww ] [ uvw vvw vww ] [ uww vww www ]

for j=1:np
    for i=1:ln
        R3_tensor_tot(:,:,1,i,j) = [F24(i,j) F28(i,j) F29(i,j);
                                    F28(i,j) F30(i,j) F34(i,j);
                                    F29(i,j) F34(i,j) F32(i,j)];

        R3_tensor_tot(:,:,2,i,j) = [F28(i,j) F30(i,j) F34(i,j);
                                    F30(i,j) F25(i,j) F31(i,j);
                                    F34(i,j) F31(i,j) F33(i,j)];

        R3_tensor_tot(:,:,3,i,j) = [F29(i,j) F34(i,j) F32(i,j);
                                    F34(i,j) F31(i,j) F33(i,j);
                                    F32(i,j) F33(i,j) F26(i,j)];
    end
end

for j=1:np
    for i=1:ln
        R3_tensor(:,:,1,i,j)=[(R3_tensor_tot(1,1,1,i,j)-3*U(i,j).*uu(i,j)-U(i,j).*U(i,j).*U(i,j)) ...
                              (R3_tensor_tot(1,2,1,i,j)-2*U(i,j).*uv(i,j)-V(i,j).*uu(i,j)-U(i,j).*U(i,j).*V(i,j)) ...
                              (R3_tensor_tot(1,3,1,i,j)-2*U(i,j).*uw(i,j)-W(i,j).*uu(i,j)-U(i,j).*U(i,j).*W(i,j));
                              (R3_tensor_tot(1,2,1,i,j)-2*U(i,j).*uv(i,j)-V(i,j).*uu(i,j)-U(i,j).*U(i,j).*V(i,j)) ...
                              (R3_tensor_tot(2,2,1,i,j)-2*V(i,j).*uv(i,j)-U(i,j).*vv(i,j)-V(i,j).*V(i,j).*U(i,j)) ...
                              (R3_tensor_tot(2,3,1,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j));
                              (R3_tensor_tot(1,3,1,i,j)-2*U(i,j).*uw(i,j)-W(i,j).*uu(i,j)-U(i,j).*U(i,j).*W(i,j)) ...
                              (R3_tensor_tot(2,3,1,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(3,3,1,i,j)-2*W(i,j).*uw(i,j)-U(i,j).*ww(i,j)-W(i,j).*W(i,j).*U(i,j))];

        R3_tensor(:,:,2,i,j)=[(R3_tensor_tot(1,1,2,i,j)-2*U(i,j).*uv(i,j)-V(i,j).*uu(i,j)-U(i,j).*U(i,j).*V(i,j)) ...
                              (R3_tensor_tot(1,2,2,i,j)-2*V(i,j).*uv(i,j)-U(i,j).*vv(i,j)-V(i,j).*V(i,j).*U(i,j)) ...
                              (R3_tensor_tot(1,3,2,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j));
                              (R3_tensor_tot(1,2,2,i,j)-2*V(i,j).*uv(i,j)-U(i,j).*vv(i,j)-V(i,j).*V(i,j).*U(i,j)) ...
                              (R3_tensor_tot(2,2,2,i,j)-3*V(i,j).*vv(i,j)-V(i,j).*V(i,j).*V(i,j)) ...
                              (R3_tensor_tot(2,3,2,i,j)-2*V(i,j).*vw(i,j)-W(i,j).*vv(i,j)-V(i,j).*V(i,j).*W(i,j));
                              (R3_tensor_tot(1,3,2,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(2,3,2,i,j)-2*V(i,j).*vw(i,j)-W(i,j).*vv(i,j)-V(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(3,3,2,i,j)-2*W(i,j).*vw(i,j)-V(i,j).*ww(i,j)-W(i,j).*W(i,j).*V(i,j))];


        R3_tensor(:,:,3,i,j)=[(R3_tensor_tot(1,1,3,i,j)-2*U(i,j).*uw(i,j)-W(i,j).*uu(i,j)-U(i,j).*U(i,j).*W(i,j)) ...
                              (R3_tensor_tot(1,2,3,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(1,3,3,i,j)-2*W(i,j).*uw(i,j)-U(i,j).*ww(i,j)-W(i,j).*W(i,j).*U(i,j));
                              (R3_tensor_tot(1,2,3,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(2,2,3,i,j)-2*V(i,j).*vw(i,j)-W(i,j).*vv(i,j)-V(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(2,3,3,i,j)-2*W(i,j).*vw(i,j)-V(i,j).*ww(i,j)-W(i,j).*W(i,j).*V(i,j));
                              (R3_tensor_tot(1,3,3,i,j)-2*W(i,j).*uw(i,j)-U(i,j).*ww(i,j)-W(i,j).*W(i,j).*U(i,j)) ...
                              (R3_tensor_tot(2,3,3,i,j)-2*W(i,j).*vw(i,j)-V(i,j).*ww(i,j)-W(i,j).*W(i,j).*V(i,j)) ...
                              (R3_tensor_tot(3,3,3,i,j)-3*W(i,j).*ww(i,j)-W(i,j).*W(i,j).*W(i,j))];
    end
end

for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        aabc(1:3,1:3,1:3)=0;
        adef=R3_tensor(:,:,:,i,j);

        for aa=1:3
        for bb=1:3
        for cc=1:3
        for dd=1:3
        for ee=1:3
        for ff=1:3
            aabc(aa,bb,cc)=aabc(aa,bb,cc)+Ralph(aa,dd)*Ralph(bb,ee) ...
                *Ralph(cc,ff)*adef(dd,ee,ff);
        end
        end
        end
        end
        end
        end

        out(j).uuu(i,1)=aabc(1,1,1);
        out(j).vvv(i,1)=aabc(2,2,2);
        out(j).www(i,1)=aabc(3,3,3);
        out(j).uuv(i,1)=aabc(1,2,1);
        out(j).uuw(i,1)=aabc(1,3,1);
        out(j).uvv(i,1)=aabc(2,2,1);
        out(j).vvw(i,1)=aabc(2,3,2);
        out(j).uww(i,1)=aabc(3,3,1);
        out(j).vww(i,1)=aabc(3,3,2);
        out(j).uvw(i,1)=aabc(2,3,1);

    end
end


%% added because of the FIK identity

%Assign second derivatives
for j=1:np
    for i=1:ln
        out(j).D45(i,1)=D45(i,j); % d2Udx2 (not rotated, mixed derivatives missing)
    end
end

%Derivatives of the Reynolds stress tensor. Tensor of Rank 3
R3_tensor(1:3,1:3,1:3,1:ln,1:np)=0;

%Form of the tensor.
% [ duudx duvdx duwdx ] [ duudy duvdy duwdy ] [ duudz duvdz duwdz ]
% [ duvdx dvvdx dvwdx ] [ duvdy dvvdy dvwdy ] [ duvdz dvvdz dvwdz ]
% [ duwdx dvwdx dwwdx ] [ duwdy dvwdy dwwdy ] [ duwdz dvwdz dwwdz ]

for j=1:np
    for i=1:ln
        R3_tensor(:,:,1,i,j) = [duudx(i,j) duvdx(i,j) duwdx(i,j);
                                duvdx(i,j) dvvdx(i,j) dvwdx(i,j);
                                duwdx(i,j) dvwdx(i,j) dwwdx(i,j)];

        R3_tensor(:,:,2,i,j) = [duudy(i,j) duvdy(i,j) duwdy(i,j);
                                duvdy(i,j) dvvdy(i,j) dvwdy(i,j);
                                duwdy(i,j) dvwdy(i,j) dwwdy(i,j)];

        R3_tensor(:,:,3,i,j) = [duudz(i,j) duvdz(i,j) duwdz(i,j);
                                duvdz(i,j) dvvdz(i,j) dvwdz(i,j);
                                duwdz(i,j) dvwdz(i,j) dwwdz(i,j)];
    end
end

for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        aabc(1:3,1:3,1:3)=0;
        adef=R3_tensor(:,:,:,i,j);

        for aa=1:3
        for bb=1:3
        for cc=1:3
        for dd=1:3
        for ee=1:3
        for ff=1:3
            aabc(aa,bb,cc)=aabc(aa,bb,cc)+Ralph(aa,dd)*Ralph(bb,ee) ...
                *Ralph(cc,ff)*adef(dd,ee,ff);
        end
        end
        end
        end
        end
        end

        out(j).duudx(i,1)=aabc(1,1,1); 
        out(j).duudy(i,1)=aabc(2,1,1); 
        out(j).dvvdx(i,1)=aabc(1,2,2); 
        out(j).dvvdy(i,1)=aabc(2,2,2); 
        out(j).duvdx(i,1)=aabc(1,2,1); 
        out(j).duvdy(i,1)=aabc(2,2,1); 
        out(j).duwdx(i,1)=aabc(1,3,1); 
        out(j).duwdy(i,1)=aabc(2,3,1); 
        out(j).dvwdx(i,1)=aabc(1,3,2); 
        out(j).dvwdy(i,1)=aabc(2,3,2); 

    end
end


%% added to test continuity
for j=1:np
    for i=1:ln
        out(j).D1(i,1)=D1(i,j); % 1. dU/dx
        out(j).D2(i,1)=D2(i,j); % 2. dU/dy
        out(j).D3(i,1)=D3(i,j); % 3. dV/dx
        out(j).D4(i,1)=D4(i,j); % 4. dV/dy
    end
end


%% YW: 20240607: Energy Flux for assessing the control gain 
%%% We take the F12, F13, F14
%% Modified Dec 09: Subtract the mean!! 
disp('YW: pu, pv, pw')
for j=1:np
    Ralph = squeeze(Ralphap(j,:,:));
    for i=1:ln
        prod=Ralph*[F12(i,j)-F1(i,j)*F4(i,j); 
                    F13(i,j)-F2(i,j)*F4(i,j); 
                    F14(i,j)-F2(i,j)*F4(i,j)];
        out(j).pu(i,1)=prod(1);
        out(j).pv(i,1)=prod(2);
        out(j).pw(i,1)=prod(3);
    end
end



end

disp('rotation done')

%%

end