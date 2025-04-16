function [wing_data] = wing_BLapprox(wing_data,istart,iend,step_spline)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nu=wing_data(1).nu;
ny=length(wing_data(1).yn);
ny
yn=wing_data(1).yn;

for iprof=istart:iend
    wing_data(iprof).P0=wing_data(iprof).P(1);
    wing_data(iprof).dPdx0=wing_data(iprof).dPdx(1);
    wing_data(iprof).DeltaPSw=zeros([1,ny]);
    wing_data(iprof).xa=wing_data(iprof).x(1);
    wing_data(iprof).ya=wing_data(iprof).y(1);
end

xa=[wing_data(istart:iend).xa];
ya=[wing_data(istart:iend).ya];
DeltaS=sqrt((xa(2:end)-xa(1:end-1)).^2+(ya(2:end)-ya(1:end-1)).^2);
nS = length(xa);
xS = zeros([1,nS]);

for iS=2:nS
    xS(iS)=xS(iS-1)+DeltaS(iS-1);
end

P0=[wing_data(istart:iend).P0];
[dPdx0S,~]=splinediff6(P0,xS,step_spline);

iS=0;
for iprof=istart:iend 
    iS=iS+1;
    wing_data(iprof).dPdx0S=dPdx0S(iS);
end

for iline=1:ny
    xaL=zeros([1,nS]);
    yaL=zeros([1,nS]);
    xSL=zeros([1,nS]);
    PP0=zeros([1,nS]);
    PL=zeros([1,nS]);
    uuL=zeros([1,nS]);
    uvL=zeros([1,nS]);
    UL=zeros([1,nS]);
    VL=zeros([1,nS]);
    dUdxL=zeros([1,nS]);
    dVdxL=zeros([1,nS]);
        
    iS=0;
    for iprof=istart:iend
        iS=iS+1;
        xaL(iS)=wing_data(iprof).x(iline);
        yaL(iS)=wing_data(iprof).y(iline);
        PP0(iS)=wing_data(iprof).P0-wing_data(iprof).P(iline);
        PL(iS)=wing_data(iprof).P(iline);
        uuL(iS)=wing_data(iprof).uu(iline);
        uvL(iS)=wing_data(iprof).uv(iline);
        UL(iS)=wing_data(iprof).U(iline);
        VL(iS)=wing_data(iprof).V(iline);
        dUdxL(iS)=wing_data(iprof).dUdx(iline);
        dVdxL(iS)=wing_data(iprof).dVdx(iline);
    end
    DeltaS=sqrt((xaL(2:end)-xaL(1:end-1)).^2+(yaL(2:end)-yaL(1:end-1)).^2);
    for iS=2:nS
        xSL(iS)=xSL(iS-1)+DeltaS(iS-1);
    end
    
    [dPdxS,~]=splinediff6(PL,xSL,step_spline);
    [dPP0dxS,~]=splinediff6(PP0,xSL,step_spline);
    [duudxS,~]=splinediff6(uuL,xSL,step_spline);
    [duvdxS,~]=splinediff6(uvL,xSL,step_spline);
    [dUdxS,~]=splinediff6(dUdxL,xSL,step_spline);
    [dVdxS,~]=splinediff6(dVdxL,xSL,step_spline);
    [d2Udx2S,~]=splinediff6(dUdxL,xSL,step_spline);
    [d2Vdx2S,~]=splinediff6(dVdxL,xSL,step_spline);
    
    iS=0;
    for iprof=istart:iend
        iS=iS+1;
        wing_data(iprof).dPdxS(iline)=dPdxS(iS);
        wing_data(iprof).DeltaPSw(iline)=dPP0dxS(iS);
        wing_data(iprof).duudxS(iline)=duudxS(iS);
        wing_data(iprof).duvdxS(iline)=duvdxS(iS);
        wing_data(iprof).dUdxS(iline)=dUdxS(iS);
        wing_data(iprof).dVdxS(iline)=dVdxS(iS);
        wing_data(iprof).d2Udx2S(iline)=d2Udx2S(iS);
        wing_data(iprof).d2Vdx2S(iline)=d2Vdx2S(iS);
    end
    
end

iy1=3;
iy2=ny;

for iprof=istart:iend
    
%     disp(iprof)
    U = zeros([ny,1]);
    V = zeros([ny,1]);
    
    dPdx = zeros([ny,1]);
    dPdy = zeros([ny,1]);
    dUdx = zeros([ny,1]);
    dUdy = zeros([ny,1]);
    dVdy = zeros([ny,1]);
    dVdx = zeros([ny,1]);
    dWdz = zeros([ny,1]);
        
    U(:)    = wing_data(iprof).U(:);
    V(:)    = wing_data(iprof).V(:);
    dPdx(:) = wing_data(iprof).dPdx(:);
    dPdy(:) = wing_data(iprof).dPdy(:);
    dUdx(:) = wing_data(iprof).dUdx(:);
    dUdy(:) = wing_data(iprof).dUdy(:);
    dVdy(:) = wing_data(iprof).dVdy(:);
    dVdx(:) = wing_data(iprof).dVdx(:);
    
    dPdxS=wing_data(iprof).dPdxS;
    d2Udx2S=wing_data(iprof).d2Udx2S;
    d2Vdx2S=wing_data(iprof).d2Vdx2S;
    duudxS=wing_data(iprof).duudxS;
    duvdxS=wing_data(iprof).duvdxS;
    
    dPdx0 = wing_data(iprof).dPdx(1);
    
    duvdy=diff6(wing_data(iprof).uv,yn);
    dvvdy=diff6(wing_data(iprof).vv,yn);
    
    d2Udy2=diff6(dUdy,yn);
    d2Vdy2=diff6(dVdy,yn);
    
    wing_data(iprof).Omegaz=(wing_data(iprof).dVdx-wing_data(iprof).dUdy) * wing_data(iprof).Ret; 

    testC1=dUdx+dVdy;
%     testNSX=U.*dUdx+V.*dUdy+dPdx-nu*d2Udy2+duvdy; % missing nu*d2Udx2 & duudx
%     testNSY=U.*dVdx+V.*dVdy+dPdy-nu*d2Vdy2+dvvdy; % missing nu*d2Vdx2 & duvdx
%     testNSX=U.*dUdx+V.*dUdy+dPdx-nu*d2Udy2+duudx+duvdy; % missing nu*d2Udx2
%     testNSY=U.*dVdx+V.*dVdy+dPdy-nu*d2Vdy2+duvdx+dvvdy; % missing nu*d2Vdx2
    testNSX=U.*dUdx+V.*dUdy+dPdx-nu*d2Udy2-nu*d2Udx2S'+duudxS'+duvdy; 
    testNSY=U.*dVdx+V.*dVdy+dPdy-nu*d2Vdy2-nu*d2Vdx2S'+duvdxS'+dvvdy; 
    testNSXS=U.*dUdx+V.*dUdy+dPdxS'-nu*d2Udy2-nu*d2Udx2S'+duudxS'+duvdy; 
    testTSL=U.*dUdx+V.*dUdy+dPdx0-nu*d2Udy2+duvdy;
    testTSLS=U.*dUdx+V.*dUdy+wing_data(iprof).dPdx0S-nu*d2Udy2+duvdy;

    wing_data(iprof).testUdUdx=U.*dUdx;
    wing_data(iprof).testVdUdy=V.*dUdy;
    wing_data(iprof).testnud2Udy2=-nu*d2Udy2;
    wing_data(iprof).testduvdy=duvdy;
    
    wing_data(iprof).testUdVdx=U.*dVdx;
    wing_data(iprof).testVdVdy=V.*dVdy;
    wing_data(iprof).testnud2Vdy2=-nu*d2Vdy2;
    wing_data(iprof).testdvvdy=dvvdy;
    
    [~,iI]=max(U); 
    iI 
    iedge=min([wing_data(iprof).i99,iI]);
    iedge=ny;
    ist=4;
    if wing_data(iprof).i99>iI 
        disp(['badness is happening, iUmax=',...
            num2str(iI),'& iedge=',num2str(iedge)])
    end 

    Uedge=U(iedge);
    % wing_data(iprof).Uedge=Uedge;
    % testC1
    
    wing_data(iprof).testC1I=trapz(yn(ist:iedge),testC1(ist:iedge));
    disp('testC1I')
    wing_data(iprof).testC1I
    wing_data(iprof).testNSXI=trapz(yn(1:iedge),testNSX(1:iedge))/(0.5*Uedge^2);
    wing_data(iprof).testNSYI=trapz(yn(1:iedge),testNSY(1:iedge))/(0.5*Uedge^2);
    wing_data(iprof).testNSXIS=trapz(yn(1:iedge),testNSXS(1:iedge))/(0.5*Uedge^2);
    wing_data(iprof).testTSLI=trapz(yn(1:iedge),testTSL(1:iedge))/(0.5*Uedge^2);
    wing_data(iprof).testTSLIS=trapz(yn(1:iedge),testTSLS(1:iedge))/(0.5*Uedge^2);
    
    wing_data(iprof).testC1=testC1;
    wing_data(iprof).testNSX=testNSX/(0.5*Uedge^2); 
    wing_data(iprof).testNSY=testNSY/(0.5*Uedge^2); 
    wing_data(iprof).testNSXS=testNSXS/(0.5*Uedge^2); 
    wing_data(iprof).testTSL=testTSL/(0.5*Uedge^2);
    wing_data(iprof).testTSLS=testTSLS/(0.5*Uedge^2);
    
end




end

