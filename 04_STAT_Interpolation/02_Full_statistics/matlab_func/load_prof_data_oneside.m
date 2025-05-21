function out = load_prof_data_oneside(side, path_mat,istart,iend,budget_flag )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% istart=30, iend=79 for low res profiles
% istart=147, iend=399 for high res profiles

disp(['load: ',path_mat])

load(path_mat);

Uinf=1.0;
rho=out.rho;
mu=out.mu;  
nu=out.nu;

np=length(out);
ln=length(out(1).U);
disp(['ln : ',num2str(ln)])
disp(['np : ',num2str(np)])

%Compute flow quantities
for i=istart:iend;

    % if i==60
    % out(i).dUdy(1) = out(i).dUdy(1)*1.05
    % end 

    out(i).tauw=mu*out(i).dUdy(1);
    out(i).dUdyw=out(i).dUdy(1);
    out(i).ut=sqrt(abs(out(i).tauw)/rho);
    out(i).Up=out(i).U/out(i).ut;
    out(i).Vp=out(i).V/out(i).ut;
    out(i).yp=out(i).yn*out(i).ut/nu;
    out(i).dUpdyp=nu/out(i).ut^2*out(i).dUdy;
end

%Compute RMSs
for i=istart:iend;
    out(i).uup=out(i).uu/out(i).ut^2;
    out(i).vvp=out(i).vv/out(i).ut^2;
    out(i).wwp=out(i).ww/out(i).ut^2;
    out(i).uvp=out(i).uv/out(i).ut^2;
end

%Use diagnostic plot to compute BL edge and integral quantities
%for 
for i=istart:iend;
% for i=147:399
    [out(i).deltap99,out(i).Uinfp,out(i).i99,out(i).H,out(i).deltasp, ...
        out(i).thetap]=diagnostic(out(i).yp',out(i).Up,sqrt(out(i).uup));
    out(i).delta99=out(i).deltap99*nu/out(i).ut;
    out(i).deltas=out(i).deltasp*nu/out(i).ut;
    out(i).dsoverd99=1/out(i).delta99*out(i).deltas;
    out(i).theta=out(i).thetap*nu/out(i).ut;
    out(i).Uinf=out(i).Uinfp*out(i).ut;
    out(i).Ret=out(i).deltap99;
    out(i).Reds=out(i).Uinf*out(i).deltas/nu;
    out(i).Reth=out(i).Uinf*out(i).theta/nu;

end

%Compute flow quantities
for i=1:length(out)
    out(i).Ue=out(i).U(out(i).i99);
    out(i).Ve=out(i).V(out(i).i99);
end

%Compute pressure gradient parameter
for i=istart:iend
% for i=147:399
    out(i).Delta=out(i).Uinfp*out(i).deltas;
    out(i).dUinfdx=out(i).dUdx(out(i).i99);
    out(i).beta=out(i).deltas/out(i).tauw*out(i).dPdx(out(i).i99);
    out(i).dPdxe=out(i).dPdx(out(i).i99);
    out(i).dtauw=out(i).deltas/out(i).tauw;
end

%compute cf:
for i=istart:iend
% for i=147:399
    switch side
        case 'pressure'
            out(i).Cf=2*(out(i).ut/out(i).Ue)^2;
        otherwise
            if out(i).dUdyw>=0
                % out(i).Cf=2*(out(i).ut/out(i).Uinf)^2;
                % out(i).Cf=2*out(i).ut/(out(i).Uinf^2);
                out(i).Cf=2*(out(i).ut/out(i).Ue)^2;    
            else
                out(i).Cf=-2*(out(i).ut/out(i).Ue)^2;
                % out(i).Cf=-2*(out(i).ut/out(i).Uinf)^2; % For global Cf
            end
    end
end

% for i = istart:iend
%     out(i).Cp = 2*(out(i).P(1)-pinf)/(rho*Uinf^2);
%     bot(i).Cp = 2*(bot(i).P(1)-pinf)/(rho*Uinf^2);
% end

for i = istart:iend
    out(i).Pw=out(i).P(1);
    out(i).Cp = 2*(out(i).P(1)-out(i).pinf)/(rho*Uinf^2);
end

%compute cf Inf:
for i = istart:iend
    switch side
        case 'pressure'
            out(i).Cfinf=2*(out(i).ut/Uinf)^2;
        otherwise
            if out(i).dUdyw>=0
                out(i).Cfinf=2*(out(i).ut/Uinf)^2;
            else
                out(i).Cfinf=-2*(out(i).ut/Uinf)^2;
            end
    end
end

%% something to better understand the FIK identity

for i = istart:iend
    
    Uinf=out(i).Uinf;
    delta99=out(i).delta99;
    
    
    eta=1-out(i).yn/delta99;
    eta(out(i).i99:end)=0;
    out(i).fiketa=eta;
    out(i).fiketa2=eta.^2;
    
    out(i).profcpf=out(i).dPdx/Uinf^2*delta99.*out(i).fiketa2;  
    
end

%% compute the Stevenson-law scaled profile
for prof = istart:iend
    u_plus=out(prof).Up;
    vwall_plus=out(prof).V(1)/out(prof).ut;
    stev = 2*u_plus./(sqrt(1+vwall_plus*u_plus)+1);
    out(prof).stev=stev;
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
        out(prof).Pk=(out(prof).Pxx+out(prof).Pyy+out(prof).Pzz);
                
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
        out(prof).Dk=(out(prof).Dxx+out(prof).Dyy+out(prof).Dzz);
        
        % MA: from the forcing terms:
        % 64. <fx*ux> % F64
        % 65. <fy*uy> % F65
        % 66. <fz*uz> % F66
        out(prof).DRTk=(out(prof).DRTx+out(prof).DRTy+out(prof).DRTz);
                
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
        out(prof).Tk=(out(prof).Txx+out(prof).Tyy+out(prof).Tzz);
        
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
        out(prof).VDk=(out(prof).VDxx+out(prof).VDyy+out(prof).VDzz);
        
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
        out(prof).VPGk=(out(prof).Pixx+out(prof).Piyy+out(prof).Pizz);
        
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
        out(prof).Ck=-(out(prof).Cxx+out(prof).Cyy+out(prof).Czz);
        
        
        out(prof).Pkp=(out(prof).Pxx+out(prof).Pyy+out(prof).Pzz)*nu/out(prof).ut^4/2;
        out(prof).Dkp=(out(prof).Dxx+out(prof).Dyy+out(prof).Dzz)*nu/out(prof).ut^4/2;
        out(prof).DRTkp=(out(prof).DRTx+out(prof).DRTy+out(prof).DRTz)*nu/out(prof).ut^4/2;
        out(prof).Tkp=(out(prof).Txx+out(prof).Tyy+out(prof).Tzz)*nu/out(prof).ut^4/2;
        out(prof).VDkp=(out(prof).VDxx+out(prof).VDyy+out(prof).VDzz)*nu/out(prof).ut^4/2;
        out(prof).VPGkp=(out(prof).Pixx+out(prof).Piyy+out(prof).Pizz)*nu/out(prof).ut^4/2;
        out(prof).Ckp=-(out(prof).Cxx+out(prof).Cyy+out(prof).Czz)*nu/out(prof).ut^4/2;
              
    end
end

%%

end 