function [outdata] = export_drag_decomposition(buffdata,istart,iend)
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
    outdata(iprof).pinf=buffdata(iprof+istart).pinf;
    
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
    
    % Edge Velocity   
    outdata(iprof).Ue=buffdata(iprof+istart).Ue;   
    
    % Shear stress and Friction Velocity 
    outdata(iprof).tauw=buffdata(iprof+istart).tauw;   
    outdata(iprof).ut=buffdata(iprof+istart).ut;   
    outdata(iprof).delta99=buffdata(iprof+istart).delta99;   
    outdata(iprof).theta=buffdata(iprof+istart).theta;   
    
    outdata(iprof).deltas=buffdata(iprof+istart).deltas;   
    outdata(iprof).Reds=buffdata(iprof+istart).Reds;   
    outdata(iprof).dsoverd99=buffdata(iprof+istart).dsoverd99;   

    
    % Clausure pressure-gradient coefficient, Reth, Re_tau, Skin-friction
    outdata(iprof).beta=buffdata(iprof+istart).beta;   
    outdata(iprof).Reth=buffdata(iprof+istart).Reth;   
    outdata(iprof).Ret=buffdata(iprof+istart).Ret;   
    outdata(iprof).Cf=buffdata(iprof+istart).Cf;   
    outdata(iprof).Cfinf=buffdata(iprof+istart).Cfinf;   
    outdata(iprof).Cp=buffdata(iprof+istart).Cp;   
    outdata(iprof).H12=buffdata(iprof+istart).deltas/buffdata(iprof+istart).theta;   

    % Inner-scaled Profile
    outdata(iprof).yp=buffdata(iprof+istart).yp;   
    outdata(iprof).Up=buffdata(iprof+istart).Up;   
    outdata(iprof).Vp=buffdata(iprof+istart).Vp;   

    % inner-scaled RMS
    outdata(iprof).uup=buffdata(iprof+istart).uup;   
    outdata(iprof).vvp=buffdata(iprof+istart).vvp;   
    outdata(iprof).wwp=buffdata(iprof+istart).wwp;   
    outdata(iprof).uvp=buffdata(iprof+istart).uvp;   

    % Turbulent Budget
    %% Production
    outdata(iprof).Pkp = buffdata(iprof+istart).Pkp;
    
    %%% Further decompose it to analyse TKE 
    outdata(iprof).Pxx = buffdata(iprof+istart).Pxx;
    outdata(iprof).Pyy = buffdata(iprof+istart).Pyy;

    %% Dissapation
    outdata(iprof).Dkp = buffdata(iprof+istart).Dkp;
    %% Dissapation Forcing term
    % outdata(iprof).DRTkp = buffdata(iprof+istart).DRTkp
    %% Turblent Transportation
    outdata(iprof).Tkp = buffdata(iprof+istart).Tkp;
    %% Viscous Diffusion 
    outdata(iprof).VDkp = buffdata(iprof+istart).VDkp;
    %% Pressure Transport
    outdata(iprof).VPGkp = buffdata(iprof+istart).VPGkp;
    %% Convection
    outdata(iprof).Ckp = buffdata(iprof+istart).Ckp;



    % Energy Gain 
    % vvv for compute the saving energy and energy-gain for opposition control  
    outdata(iprof).vvv=buffdata(iprof+istart).vvv; 
    % Energy Flux for computing the energy gain for opposition control 
    outdata(iprof).pu=buffdata(iprof+istart).pu; 
    outdata(iprof).pv=buffdata(iprof+istart).pv; 
    outdata(iprof).pw=buffdata(iprof+istart).pw; 

    % Velocity profiles by Stevenson Law 
    outdata(iprof).stev=buffdata(iprof+istart).stev; 

     %% Aerodynamic Efficiency 
    outdata(iprof).Cl=buffdata(iprof+istart).cl2;
    outdata(iprof).Cd=buffdata(iprof+istart).cd2;
    outdata(iprof).Cd_p=buffdata(iprof+istart).cd_p;
    outdata(iprof).Cd_f=buffdata(iprof+istart).cd_f;
    outdata(iprof).LD=buffdata(iprof+istart).LD;

end

end

