function [outdata] = export_drag_decomposition_wake(buffdata,istart,iend)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nprof=iend-istart;

for iprof=1:nprof;
    
    % Corrdinates and profiles
    outdata(iprof).xa=buffdata(iprof+istart).xa;
    outdata(iprof).ya=buffdata(iprof+istart).ya;
    outdata(iprof).alpha=buffdata(iprof+istart).alpha;
    outdata(iprof).x=buffdata(iprof+istart).x;
    outdata(iprof).y=buffdata(iprof+istart).y;
    outdata(iprof).yn=buffdata(iprof+istart).yn;
    % Characteristics
    outdata(iprof).mu=buffdata(iprof+istart).mu;
    outdata(iprof).nu=buffdata(iprof+istart).nu;
    outdata(iprof).rho=buffdata(iprof+istart).rho;
    
    % Velocity
    outdata(iprof).U=buffdata(iprof+istart).U;
    outdata(iprof).V=buffdata(iprof+istart).V;
    outdata(iprof).W=buffdata(iprof+istart).W;
    
    %Pressure and pressure fluctuations 
    outdata(iprof).P=buffdata(iprof+istart).P; 
    outdata(iprof).pp=buffdata(iprof+istart).pp; 
    
    % Pressure gradient
    outdata(iprof).dPdx=buffdata(iprof+istart).dPdx;
    outdata(iprof).dPdy=buffdata(iprof+istart).dPdy;
    outdata(iprof).dPdz=buffdata(iprof+istart).dPdz;
    
    % Fluctuation/Stress 
    outdata(iprof).uu=buffdata(iprof+istart).uu;
    outdata(iprof).vv=buffdata(iprof+istart).vv;
    outdata(iprof).ww=buffdata(iprof+istart).ww;
    outdata(iprof).uv=buffdata(iprof+istart).uv;
    outdata(iprof).uw=buffdata(iprof+istart).uw;
    outdata(iprof).vw=buffdata(iprof+istart).vw;
    
    % Velocity Gradient 
    outdata(iprof).dUdx=buffdata(iprof+istart).dUdx;
    outdata(iprof).dUdy=buffdata(iprof+istart).dUdy;
    outdata(iprof).dVdx=buffdata(iprof+istart).dVdx;    
    outdata(iprof).dVdy=buffdata(iprof+istart).dVdy;   
    
    % Freestream Velocity 
    outdata(iprof).Uinf=buffdata(iprof+istart).Uinf;
    
    % Turbulent Budget
    %% Production
    outdata(iprof).Pk = buffdata(iprof+istart).Pk;
    
    %%% Further decompose it to analyse TKE 
    outdata(iprof).Pxx = buffdata(iprof+istart).Pxx;
    outdata(iprof).Pyy = buffdata(iprof+istart).Pyy;

    %% Dissapation
    outdata(iprof).Dk = buffdata(iprof+istart).Dk;
    %% Dissapation Forcing term
    % outdata(iprof).DRTk = buffdata(iprof+istart).DRTk
    %% Turblent Transportation
    outdata(iprof).Tk = buffdata(iprof+istart).Tk;
    %% Viscous Diffusion 
    outdata(iprof).VDk = buffdata(iprof+istart).VDk;
    %% Pressure Transport
    outdata(iprof).VPGk = buffdata(iprof+istart).VPGk;
    %% Convection
    outdata(iprof).Ck = buffdata(iprof+istart).Ck;

    % Energy Gain 
    % vvv for compute the saving energy and energy-gain for opposition control  
    outdata(iprof).vvv=buffdata(iprof+istart).vvv; 
    % Energy Flux for computing the energy gain for opposition control 
    outdata(iprof).pu=buffdata(iprof+istart).pu; 
    outdata(iprof).pv=buffdata(iprof+istart).pv; 
    outdata(iprof).pw=buffdata(iprof+istart).pw; 

end

end

