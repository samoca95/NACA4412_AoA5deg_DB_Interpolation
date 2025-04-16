function [x_pts,y_pts,npoints,yn,alphau,alphal] = create_profiles_wing2( m,p,xx,nle,tol,dy1,yn,r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

disp('compute wing profiles')

[xu,yu,alphau,xl,yl,alphal] = naca_prof4( m,p,xx,nle,tol );

y1=0; %First point is at the wall

%Find the progressive ratio and the number of points from the wall to
%delta1. We round towards the closest higher integer.
n=ceil(log(-(y1*(r-1)+yn*(1-r)-dy1)/dy1)/log(r)+1);

%Generate sequence of points up to delta1
yn(1)=y1;
for i=2:n
    yn(i)=yn(1)+dy1*(1-r^(i-1))/(1-r);
end

[xprofu,yprofu,xprofl,yprofl]=rot(xu,yu,xl,yl,yn,alphau,alphal);

npoints=0;

%Rearrange interpolating mesh
%Upper surface
for i=1:length(xprofu(1,:))
    for j=1:length(xprofu(:,1))
        npoints=npoints+1;
        x_pts(npoints,1)=xprofu(j,i);
        y_pts(npoints,1)=yprofu(j,i);
    end
end

%Lower surface
for i=2:length(xprofl(1,:))
    for j=1:length(xprofl(:,1))
        npoints=npoints+1;
        x_pts(npoints,1)=xprofl(j,i);
        y_pts(npoints,1)=yprofl(j,i);
    end
end

%Store wall-normal profile yn
for i=1:length(yn)
    npoints=npoints+1;
    x_pts(npoints,1)=-1;
    y_pts(npoints,1)=yn(i);
end
end

