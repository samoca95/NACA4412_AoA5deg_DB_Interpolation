function out = load_prof_data_wake( path_mat,istart,iend,budget_flag )
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

%% budget turbulent-kinetic energy 
% 
if budget_flag==1
    for prof=istart:iend

        out(prof).Uinf=Uinf;

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
        
    end
end

end 