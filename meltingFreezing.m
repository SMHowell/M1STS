%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M,OUT] = meltingFreezing(M,IN,COMP,OUT)


%%%%%%%%%%%%%%%%%%
% Ocean Melting/Freezing
%%%%%%%%%%%%%%%%%%

% Calculate heat extracted from ocean
dE_ocn = (M.T(M.iOcnTop)-M.Tm_ocn)*M.rho(M.iOcnTop)*M.Cp(M.iOcnTop)*(4/3)*pi*(M.r_s(M.iOcnTop)^3 - M.r_s(M.iOcnTop-1)^3);

% Also, check whether the interface has warmed above Tm
dE_ice = (M.T(M.iOcnTop+1)-M.Tm_ocn)*M.rho(M.iOcnTop+1)*M.Cp(M.iOcnTop+1)*(4/3)*pi*(M.r_s(M.iOcnTop+1)^3 - M.r_s(M.iOcnTop)^3);
dE_ice(dE_ice<0) = 0;         % Ice is allowed to be colder than melting
if M.T(M.iOcnTop+1)>M.Tm_ocn % Ice is not allowed to be warmer
    M.T(M.iOcnTop+1) = M.Tm_ocn;
end

% Calculate change in melt
dE = dE_ocn + dE_ice;
dm = dE/IN.L;    % Total change in melt mass
dv = dm/M.rhoOcn; % Total change in melt volume

M.dm_dt = dm/M.dt;

% New interface location
M.rOcnTop = ((3/(4*pi))*dv+M.rOcnTop^3)^(1/3);
M.iOcnTop = find((M.rOcnTop - M.r>0)>0,1,'last');

% Set melt fractions
M.vfm(M.r_s<M.rOcnTop) = 1;
M.vfm(M.r_s>M.rOcnTop) = 0;
M.vfm(M.iOcnTop)       = (4/3)*pi*(M.rOcnTop^3 - M.r(M.iOcnTop)^3)/M.V_s(M.iOcnTop);

% Restore correct temperatures now that energy is back in the right place
M.T(M.T >M.Tm_ocn) = M.Tm_ocn;
M.T(M.r<=M.rOcnTop) = M.Tm_ocn;

% Calculate change in melt
dE = dE_ocn + dE_ice;
dm = dE/IN.L;    % Total change in melt mass
dv = dm/M.rhoOcn; % Total change in melt volume

M.dm_dt = dm/M.dt;

% New interface location
M.rOcnTop = ((3/(4*pi))*dv+M.rOcnTop^3)^(1/3);
M.iOcnTop = find((M.rOcnTop - M.r>0)>0,1,'last');

% Set melt fractions
M.vfm(M.r_s<M.rOcnTop) = 1;
M.vfm(M.r_s>M.rOcnTop) = 0;
M.vfm(M.iOcnTop) = (4/3)*pi*(M.rOcnTop^3 - M.r(M.iOcnTop)^3)/M.V_s(M.iOcnTop);

% Restore correct temperatures now that energy is back in the right place
M.T(M.r<=M.rOcnTop) = M.Tm_ocn;


%%%%%%%%%%%%%%%%%%
% Reservoir Melting/Freezing
%%%%%%%%%%%%%%%%%%
% NOTE: This works by getting the energy extracted from an annulus of water
% spanning the globe. We only use this energy value to calculate change in
% reservoir interface position. By calculating change in reservoir radius
% at the top and bottom, we can estimate overall change in reservoir volume
% effectively.

if M.vRes>0  

    % Grab reservoir properties from composition data
    [T, rho, vlf, vif] = coolFreeze(M, COMP, IN);

    % Update parameters
    M.vif       = vif;
    M.rhoRes    = rho;
    M.fV        = 1 - vlf;
    M.Tm_res    = T;
    OUT.Tmelt   = T;
    
    % Change in reservoir size
    vResOld     = M.vRes;
    M.vRes      = (4/3)*pi*IN.rRes^3 * vlf;     % New reservoir volume
    vResDiff    = M.vRes - vResOld;
    vDiffTop    = vResDiff * M.ratioEtop;       % Volume diff at res top
    vDiffBot    = vResDiff * M.ratioEbot;       % Volume diff at res bottom (not sure these are very clean...)
        
    % Change in interface location
    M.rResTop = ((3/(4*pi))*vDiffTop+M.rResTop^3)^(1/3);
    M.rResBot = (-(3/(4*pi))*vDiffBot+M.rResBot^3)^(1/3);

    % Freezing front velocity (average top and bottom)
    M.freezeRate = ((vResOld*3/4/pi)^(1/3)-(M.vRes*3/4/pi)^(1/3)) / M.dt * 100/365.25/24/3600; % in cm/yr
     
    % Change in reservoir size and frozen fraction
    M.rRes  = (M.rResTop-M.rResBot)/2; % Reservoir radius
    
    M.iResTop = find((M.rResTop - M.r>0)>0,1,'last'); % Reservoir top interface element index
    M.iResBot = find((M.rResBot - M.r>0)>0,1,'last'); % Reservoir bottom interface element index

    % Set Temp
    M.T(M.iResBot:M.iResTop) = M.Tm_res;
    
end
end
