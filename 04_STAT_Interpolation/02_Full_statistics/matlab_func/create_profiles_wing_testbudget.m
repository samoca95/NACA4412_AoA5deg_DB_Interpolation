function [x_pts,y_pts,npoints,yn,alphau,alphal,x] = create_profiles_wing_testbudget( m,p,xx,dx,nle,tol,dy1,yn,ry)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

disp('compute wing profiles')

c=1;       %Chord length

%Between x/c=0 and 0.3 we need a progressive ratio, clustering many points
%close to the leading edge, and progressively converging towards the
%uniform value dx.


%We fix the first point as 0, the last one as 0.3, and the spacing of the
%last point is dx. Knowing the number of points, we can define a geometric
%series to solve for the spacing of the first point and the progressive
%ratio
x1=0.0; %Location of first point (leading edge)
xn=0.295; %Location of last point (after which uniform dx is considered)
xn1=xn-dx; %Last spacing is uniform value of dx



%From sum of power series, it is possible to solve for the progression
%ratio between spacings r
find_r=@(r) r^(1-nle)*(xn-xn1)*(r^nle-r)/(r-1)-xn+x1;
r=fzero(find_r,1.15);

%Then we solve for the spacing between the first and second points
dx1=r^(2-nle)*(xn-xn1);

%Generate sequence of points for leading edge
x(1)=x1;
for i=2:nle
    x(i)=x(1)+dx1*(1-r^(i-1))/(1-r);
end

%Generate sequence of points for full wing
x=[x(1:end-1) x(end):dx:c];

[xu,yu,alphau,xl,yl,alphal] = naca_prof2( m,p,xx,x,tol );

%%

y1=0; %First point is at the wall

%Find the progressive ratio and the number of points from the wall to
%delta1. We round towards the closest higher integer.
n=ceil(log(-(y1*(ry-1)+yn*(1-ry)-dy1)/dy1)/log(ry)+1);

%Generate sequence of points up to delta1
yn(1)=y1;
for i=2:n
    yn(i)=yn(1)+dy1*(1-ry^(i-1))/(1-ry);
end

[xprofu,yprofu,xprofl,yprofl]=rot(xu,yu,xl,yl,yn,alphau,alphal);

%%

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

