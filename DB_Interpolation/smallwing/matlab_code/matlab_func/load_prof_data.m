function [ top, bot ] = load_prof_data( path_mat,istart,iend,budget_flag )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% istart=30, iend=79 for low res profiles
% istart=147, iend=399 for high res profiles

disp(['load: ',path_mat])

load(path_mat);

Uinf=1.0;
rho=top.rho;
mu=top.mu;
nu=top.nu;

bot=bottom;

np=length(top);
ln=length(top(1).U);

%Compute flow quantities
for i=istart:iend
    top(i).tauw=mu*top(i).dUdy(1);
    top(i).dUdyw=top(i).dUdy(1);
    top(i).ut=sqrt(abs(top(i).tauw)/rho);
    top(i).Up=top(i).U/top(i).ut;
    top(i).Vp=top(i).V/top(i).ut;
    top(i).yp=top(i).yn*top(i).ut/nu;
    top(i).dUpdyp=nu/top(i).ut^2*top(i).dUdy;
end

for i=istart:iend
    bot(i).tauw=mu*bot(i).dUdy(1);
    bot(i).dUdyw=bot(i).dUdy(1);
    bot(i).ut=sqrt(abs(bot(i).tauw)/rho);
    bot(i).Up=bot(i).U/bot(i).ut;
    bot(i).Vp=bot(i).V/bot(i).ut;
    bot(i).yp=bot(i).yn*bot(i).ut/nu;
    bot(i).dUpdyp=nu/bot(i).ut^2*bot(i).dUdy;
end

%Compute RMSs
for i=istart:iend
    top(i).uup=top(i).uu/top(i).ut^2;
    top(i).vvp=top(i).vv/top(i).ut^2;
    top(i).wwp=top(i).ww/top(i).ut^2;
    top(i).uvp=top(i).uv/top(i).ut^2;
end

for i=istart:iend
    bot(i).uup=bot(i).uu/bot(i).ut^2;
    bot(i).vvp=bot(i).vv/bot(i).ut^2;
    bot(i).wwp=bot(i).ww/bot(i).ut^2;
    bot(i).uvp=bot(i).uv/bot(i).ut^2;
end

%Use diagnostic plot to compute BL edge and integral quantities
%for 
for i=istart:iend
% for i=147:399
    [top(i).deltap99,top(i).Uinfp,top(i).i99,top(i).H,top(i).deltasp, ...
        top(i).thetap]=diagnostic(top(i).yp,top(i).Up,sqrt(top(i).uup));
    top(i).delta99=top(i).deltap99*nu/top(i).ut;
    top(i).deltas=top(i).deltasp*nu/top(i).ut;
    top(i).dsoverd99=1/top(i).delta99*top(i).deltas;
    top(i).theta=top(i).thetap*nu/top(i).ut;
    top(i).Uinf=top(i).Uinfp*top(i).ut;
    top(i).Ret=top(i).deltap99;
    top(i).Reds=top(i).Uinf*top(i).deltas/nu;
    top(i).Reth=top(i).Uinf*top(i).theta/nu;

end

for i=istart:iend
% for i=147:399
    [bot(i).deltap99,bot(i).Uinfp,bot(i).i99,bot(i).H,bot(i).deltasp, ...
        bot(i).thetap]=diagnostic(bot(i).yp,bot(i).Up,sqrt(bot(i).uup));
    bot(i).delta99=bot(i).deltap99*nu/bot(i).ut;
    bot(i).deltas=bot(i).deltasp*nu/bot(i).ut;
    bot(i).dsoverd99=1/bot(i).delta99*bot(i).deltas;
    bot(i).theta=bot(i).thetap*nu/bot(i).ut;
    bot(i).Uinf=bot(i).Uinfp*bot(i).ut;
    bot(i).Ret=bot(i).deltap99;
    bot(i).Reds=bot(i).Uinf*bot(i).deltas/nu;
    bot(i).Reth=bot(i).Uinf*bot(i).theta/nu;

end


%Compute flow quantities
for i=1:length(top)
    top(i).Ue=top(i).U(top(i).i99);
    top(i).Ve=top(i).V(top(i).i99);
end

for i=1:length(bot)
    bot(i).Ue=bot(i).U(bot(i).i99);
    bot(i).Ve=bot(i).V(bot(i).i99);
end

%Compute pressure gradient parameter
for i=istart:iend
% for i=147:399
    top(i).Delta=top(i).Uinfp*top(i).deltas;
    top(i).dUinfdx=top(i).dUdx(top(i).i99);
    top(i).beta=top(i).deltas/top(i).tauw*top(i).dPdx(top(i).i99);
    top(i).dPdxe=top(i).dPdx(top(i).i99);
    top(i).dtauw=top(i).deltas/top(i).tauw;
end

for i=istart:iend
% for i=147:399
    bot(i).Delta=bot(i).Uinfp*bot(i).deltas;
    bot(i).dUinfdx=bot(i).dUdx(bot(i).i99);
    bot(i).beta=bot(i).deltas/bot(i).tauw*bot(i).dPdx(bot(i).i99);
    bot(i).dPdxe=bot(i).dPdx(bot(i).i99);
    bot(i).dtauw=bot(i).deltas/bot(i).tauw;
end

%compute cf:
for i=istart:iend
% for i=147:399
    if top(i).dUdyw>=0
        top(i).Cf=2*(top(i).ut/top(i).Ue)^2;
    else
        top(i).Cf=-2*(top(i).ut/top(i).Ue)^2;
    end
end

for i=istart:iend
% for i=147:399
    bot(i).Cf=2*(bot(i).ut/bot(i).Ue)^2;
end


% for i = istart:iend
%     top(i).Cp = 2*(top(i).P(1)-pinf)/(rho*Uinf^2);
%     bot(i).Cp = 2*(bot(i).P(1)-pinf)/(rho*Uinf^2);
% end

for i = istart:iend
    top(i).Pw=top(i).P(1);
    top(i).Cp = 2*(top(i).P(1)-top(i).pinf)/(rho*Uinf^2);
    bot(i).Pw=bot(i).P(1);
    bot(i).Cp = 2*(bot(i).P(1)-top(i).pinf)/(rho*Uinf^2);
end

%compute cf Inf:
for i = istart:iend
    if top(i).dUdyw>=0
        top(i).Cfinf=2*(top(i).ut/Uinf)^2;
    else
        top(i).Cfinf=-2*(top(i).ut/Uinf)^2;
    end
end

for i = istart:iend
    bot(i).Cfinf=2*(bot(i).ut/Uinf)^2;
end

%% something to better understand the FIK identity

for i = istart:iend
    
    Uinf=top(i).Uinf;
    delta99=top(i).delta99;
    
    
    eta=1-top(i).yn/delta99;
    eta(top(i).i99:end)=0;
    top(i).fiketa=eta;
    top(i).fiketa2=eta.^2;
    
    top(i).profcpf=top(i).dPdx/Uinf^2*delta99.*top(i).fiketa2;
    
    Uinf=top(i).Uinf;
    delta99=top(i).delta99;    
    eta=1-bot(i).yn/delta99;
    eta(bot(i).i99:end)=0;
    bot(i).fiketa=eta; 
    bot(i).fiketa2=eta.^2; 
    
    bot(i).profcpf=bot(i).dPdx/Uinf^2*delta99.*bot(i).fiketa2;
    
    
end

%% compute the Stevenson-law scaled profile


for prof = istart:iend
    u_plus=top(prof).Up;
    vwall_plus=top(prof).V(1)/top(prof).ut;
    stev = 2*u_plus./(sqrt(1+vwall_plus*u_plus)+1);
    top(prof).stev=stev;
end

for prof = istart:iend
    u_plus=bot(prof).Up;
    vwall_plus=bot(prof).V(1)/bot(prof).ut;
    stev = 2*u_plus./(sqrt(1+vwall_plus*u_plus)+1);
    bot(prof).stev=stev;
end

%% budget turbulent-kinetic energy 
% 
if budget_flag==1
    for prof=istart:iend

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
        top(prof).Pk=(top(prof).Pxx+top(prof).Pyy+top(prof).Pzz);
                
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
        top(prof).Dk=(top(prof).Dxx+top(prof).Dyy+top(prof).Dzz);
        
        % MA: from the forcing terms:
        % 64. <fx*ux> % F64
        % 65. <fy*uy> % F65
        % 66. <fz*uz> % F66
        top(prof).DRTk=(top(prof).DRTx+top(prof).DRTy+top(prof).DRTz);
                
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
        top(prof).Tk=(top(prof).Txx+top(prof).Tyy+top(prof).Tzz);
        
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
        top(prof).VDk=(top(prof).VDxx+top(prof).VDyy+top(prof).VDzz);
        
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
        top(prof).VPGk=(top(prof).Pixx+top(prof).Piyy+top(prof).Pizz);
        
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
        top(prof).Ck=-(top(prof).Cxx+top(prof).Cyy+top(prof).Czz);
        
        
        top(prof).Pkp=(top(prof).Pxx+top(prof).Pyy+top(prof).Pzz)*nu/top(prof).ut^4/2;
        top(prof).Dkp=(top(prof).Dxx+top(prof).Dyy+top(prof).Dzz)*nu/top(prof).ut^4/2;
        top(prof).DRTkp=(top(prof).DRTx+top(prof).DRTy+top(prof).DRTz)*nu/top(prof).ut^4/2;
        top(prof).Tkp=(top(prof).Txx+top(prof).Tyy+top(prof).Tzz)*nu/top(prof).ut^4/2;
        top(prof).VDkp=(top(prof).VDxx+top(prof).VDyy+top(prof).VDzz)*nu/top(prof).ut^4/2;
        top(prof).VPGkp=(top(prof).Pixx+top(prof).Piyy+top(prof).Pizz)*nu/top(prof).ut^4/2;
        top(prof).Ckp=-(top(prof).Cxx+top(prof).Cyy+top(prof).Czz)*nu/top(prof).ut^4/2;
        
        bot(prof).Pk=(bot(prof).Pxx+bot(prof).Pyy+bot(prof).Pzz);
        bot(prof).Dk=(bot(prof).Dxx+bot(prof).Dyy+bot(prof).Dzz);
        bot(prof).DRTk=(bot(prof).DRTx+bot(prof).DRTy+bot(prof).DRTz);
        bot(prof).Tk=(bot(prof).Txx+bot(prof).Tyy+bot(prof).Tzz);
        bot(prof).VDk=(bot(prof).VDxx+bot(prof).VDyy+bot(prof).VDzz);
        bot(prof).VPGk=(bot(prof).Pixx+bot(prof).Piyy+bot(prof).Pizz);
        bot(prof).Ck=-(bot(prof).Cxx+bot(prof).Cyy+bot(prof).Czz);
        
        bot(prof).Pkp=(bot(prof).Pxx+bot(prof).Pyy+bot(prof).Pzz)*nu/bot(prof).ut^4/2;
        bot(prof).Dkp=(bot(prof).Dxx+bot(prof).Dyy+bot(prof).Dzz)*nu/bot(prof).ut^4/2;
        bot(prof).DRTk=(bot(prof).DRTx+bot(prof).DRTy+bot(prof).DRTz)*nu/bot(prof).ut^4/2;
        bot(prof).Tkp=(bot(prof).Txx+bot(prof).Tyy+bot(prof).Tzz)*nu/bot(prof).ut^4/2;
        bot(prof).VDkp=(bot(prof).VDxx+bot(prof).VDyy+bot(prof).VDzz)*nu/bot(prof).ut^4/2;
        bot(prof).VPGkp=(bot(prof).Pixx+bot(prof).Piyy+bot(prof).Pizz)*nu/bot(prof).ut^4/2;
        bot(prof).Ckp=-(bot(prof).Cxx+bot(prof).Cyy+bot(prof).Czz)*nu/bot(prof).ut^4/2;
        
    end
end

%%




end 