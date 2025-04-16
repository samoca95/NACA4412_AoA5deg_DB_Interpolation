function [ oxu,oyu,oalphau,oxl,oyl,oalphal] = naca_prof4( m,p,xx,npoint,tol )

%% It gives the (xu,yu) and (xl,yl) points corresponding to inputx points along the chord, 
%% along with the relative angle normal to the surface, alphau and alphal.
%% + REDISTRIBUTION OF POINTS ALONG THE SUCTION AND THE PRESSURE SIDE

tic


M=m/100;
P=p/10;
T=xx/100;

beta=0:tol:pi;

x = 1/2*(1-cos(beta));

%%

% Front (0<=x<p)

xf=x(x<P);

yc_f = M/(P^2)*(2*P*xf-xf.^2);

dycdx_f = 2*M/(P^2)*(P-xf);

theta_f = atan(dycdx_f);

% Back (p<=x<=1)

xb=x(x>=P);

yc_b = M/(1-P)^2*(1-2*P+2*P*xb-xb.^2);

dycdx_b = 2*M/(1-P)^2*(P-xb);

theta_b = atan(dycdx_b);

% NACA coefficients:

a0=0.2969;
a1=-0.126;
a2=-0.3516;
a3=0.2843;
a4=-0.1036;

% thickness distribution:

yt = T/0.2*(a0*x.^(0.5)+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);

% chord and theta distribution:
xc = [xf,xb];
yc = [yc_f, yc_b];

theta=[theta_f,theta_b];

% point location:

xu = xc - yt.*sin(theta);
yu = yc + yt.*cos(theta); 

xl= xc + yt.*sin(theta);
yl= yc - yt.*cos(theta);

%%
% 
% figure
% ref=importdata('naca4412.txt');
% plot(xu,yu,'-og'); hold on
% plot(xl,yl,'-ob'); hold on
% plot(ref(:,1),ref(:,2),'-ok')
% axis equal

%%

xtmp1 = [xu, xu(1)];
ytmp1 = [yu, yu(1)];
xtmp2 = [xu(end), xu];
ytmp2 = [yu(end), yu];
nxu =   (diff(ytmp1)+diff(ytmp2))/2.0;
nyu =  -(diff(xtmp1)+diff(xtmp2))/2.0;
alphau = atan(-nxu./nyu);

xtmp1 = [xl, xl(1)];
ytmp1 = [yl, yl(1)];
xtmp2 = [xl(end), xl];
ytmp2 = [yl(end), yl];
nxl =   (diff(ytmp1)+diff(ytmp2))/2.0;
nyl =  -(diff(xtmp1)+diff(xtmp2))/2.0;
alphal = atan(-nxl./nyl);

alphau(end)=alphau(end-1);

%%

beta=0:pi/npoint:pi;
inputx = 1/2*(1-cos(beta));

[~,ixmin]=min(xu);
xL=x(ixmin);

true_suction_side=inputx>xL;
inputxu1 = inputx(true_suction_side);

inputxu2 = inputx(inputx<xL);

n_false_ss = npoint-sum(true_suction_side);

beta=0:pi/(npoint-2*(n_false_ss+1)):pi;
inputxl = 1/2*(1-cos(beta));

%%

oxu=interp1(x,xu,inputxu1,'spline');
oyu=interp1(x,yu,inputxu1,'spline');
oalphau=interp1(x,alphau,inputxu1,'spline');

false_oxu=interp1(x,xu,inputxu2,'spline');
false_oyu=interp1(x,yu,inputxu2,'spline');
false_oalphau=interp1(x,alphau,inputxu2,'spline');

temp_oxl=interp1(x,xl,inputxl,'spline');
temp_oyl=interp1(x,yl,inputxl,'spline');
temp_oalphal=interp1(x,alphal,inputxl,'spline');
temp_oalphal(1)=temp_oalphal(2);

oxl=[false_oxu,temp_oxl];
oyl=[false_oyu,temp_oyl];
oalphal=[false_oalphau,temp_oalphal];

toc

disp('done compute prof')

end


