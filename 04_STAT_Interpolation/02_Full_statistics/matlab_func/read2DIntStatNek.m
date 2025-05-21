function [ductstat] = read2DIntStatNek(pathintstat,size_xx,size_yy)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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


% Read interpolated field: int_fld file
fname         = [pathintstat,'/ZSTAT/int_fld'];
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

% Read interpolating mesh: x.fort and y.fort
fname_x         = [pathintstat,'/ZSTAT/x.fort'];
[fid_x,message_x] = fopen(fname_x,'r','ieee-le');
hdr_x             = fread(fid_x,1,'int32');
dum1              = fread(fid_x,3,'int32');
x_pts             = fread(fid_x,npoints,'*float64');
fclose(fid_x);

fname_y         = [pathintstat,'/ZSTAT/y.fort'];
[fid_y,message_y] = fopen(fname_y,'r','ieee-le');
hdr_y             = fread(fid_y,1,'int32');
dum2              = fread(fid_y,3,'int32');
y_pts             = fread(fid_y,npoints,'*float64');
fclose(fid_y);

%%

ductstat.X=reshape(x_pts,size_xx,size_yy);
ductstat.Y=reshape(y_pts,size_xx,size_yy)+1;


%Define data arrays
Mat=zeros(size_yy,size_xx);
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

    %     % Arrange current field in array form
    %     for i=1:2*np
    %         Mat(:,i)=MatR(ln*(i-1)+1:ln*(i-1)+ln);
    %     end
    Mat=MatR;

    %     if ii<=nstat
    %         % Generata variable name for field X as FX
    %         v=genvarname('F', who);
    %     else
    %         % Generata variable name for field X as DX
    %         v=genvarname('D', who);
    %     end
    %
    %     % Store matrix in variable FX or Dx
    %     evalc([v '=INOUT(Mat)']);


%% mean velocity

    if ii==1
        ductstat.U=reshape(Mat,size_xx,size_yy);
    end
    if ii==2
        ductstat.V=reshape(Mat,size_xx,size_yy);
    end
    if ii==3
        ductstat.W=reshape(Mat,size_xx,size_yy);
    end

%% velocity fluctuations

    if ii==5
        ductstat.uu=reshape(Mat,size_xx,size_yy);
    end
    if ii==6
        ductstat.vv=reshape(Mat,size_xx,size_yy);
    end
    if ii==7
        ductstat.ww=reshape(Mat,size_xx,size_yy);
    end
    if ii==9
        ductstat.uv=reshape(Mat,size_xx,size_yy);
    end
    if ii==11
        ductstat.uw=reshape(Mat,size_xx,size_yy);
    end
    if ii==10
        ductstat.vw=reshape(Mat,size_xx,size_yy);
    end

%% velocity derivatives 

    if ii==nstat+1
        ductstat.dUdx=reshape(Mat,size_xx,size_yy);
    end
    if ii==nstat+2
        ductstat.dUdy=reshape(Mat,size_xx,size_yy);
    end
    if ii==nstat+3
        ductstat.dVdx=reshape(Mat,size_xx,size_yy);
    end
    if ii==nstat+4
        ductstat.dVdy=reshape(Mat,size_xx,size_yy);
    end
    if ii==nstat+5
        ductstat.dWdx=reshape(Mat,size_xx,size_yy);
    end
    if ii==nstat+6
        ductstat.dWdy=reshape(Mat,size_xx,size_yy);
    end
%     if ii==nstat+7
%         ductstat.dPdx=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+8
%         ductstat.dPdy=reshape(Mat,size_xx,size_yy);
%     end

%% fluctuations derivatives

%     if ii==nstat+9
%         ductstat.duudx=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+10
%         ductstat.duudy=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+11
%         ductstat.dvvdx=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+12
%         ductstat.dvvdy=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+13
%         ductstat.dwwdx=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+14
%         ductstat.dwwdy=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+15
%         ductstat.dppdx=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+16
%         ductstat.dppdy=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+17
%         ductstat.duvdx=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+18
%         ductstat.duvdy=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+19
%         ductstat.dvwdx=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+20
%         ductstat.dvwdy=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+21
%         ductstat.duwdx=reshape(Mat,size_xx,size_yy);
%     end
%     if ii==nstat+22
%         ductstat.duwdy=reshape(Mat,size_xx,size_yy);
%     end


end

fclose(fid);



end