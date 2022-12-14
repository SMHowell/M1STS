%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = differentiate(M,IN,BOD)

%%%%%%%%%%%%%%%%%%
% Extract water
%%%%%%%%%%%%%%%%%%
% We have to remove the ice fraction in the rocky interior towards the
% surface. Once it reaches the melting temperature, the volume erupts,
% either adding to an existing ocean or refreezing on the seafloor and
% adding energy to the seafloor temperature.

%%%%%%%
% Set up problem
%%%%%%%
% Indices where icy differentiation should occur
[fmDiff_H2O] = s2nVolumetric(M,M.fmDiff_H2O);
ind = find(M.fV_H2O>0 & M.fV_sil>0 & M.T>IN.Tm_ocn & (1-abs(fmDiff_H2O)>1e-5));

% Energy available
dE  = 2*pi.*M.rho(ind).*M.Cp(ind).*M.r(ind).^2.*...
    (M.dr(ind-1)+M.dr(ind)).*(M.T(ind)-IN.Tm_ocn);

% Total energy required to melt water present
% Get mass fractions of differentiated material on nodes
dE_max = M.fm_sil(ind).* (1 - fmDiff_H2O(ind)).* IN.fm0_H2O.*M.V(ind).*M.rho(ind)*BOD.LH2O;

%%%%%%%
% Calculate mass change
%%%%%%%
dExcess = dE-dE_max;
dExcess(dExcess<0) = 0;

% Melted mass of water on node
dm = min(dE,dE_max)./BOD.LH2O;

if any(abs(dm)>0)
    % Because we withdrew the extra energy, make sure to set
    % temperature accordingly.
    M.T(ind) = IN.Tm_ocn;
    
    dT = dExcess./(2*pi*M.rho(ind).*M.Cp(ind).*M.r(ind).^2.*...
        (M.dr(ind-1)+M.dr(ind)));
    M.T(ind) = M.T(ind) + dT;
    
    
    %%%%%%%
    % Account for mass loss
    %%%%%%%
    % Change in rock volume on node
    dV   = -dm./BOD.rhoIce_0;
    dV_n = zeros(1,M.Nz); dV_n(ind) = dV;
    
    % Change in liquid water volume
    dV_n(M.iOcnTop) = -sum(dV_n); %.*BOD.rhoIce_0/BOD.rhoOcn;
    
    % Propogate change
    rTemp  = M.r;
    for i=ind:M.Nz
        rTemp(i) = ((3/(4*pi))*(dV_n(i)+M.V_s(i-1)) + rTemp(i-1)^3).^(1/3);
    end
    
    % Temperature update using upwind finite differences
    dr = rTemp-M.r;
    dT = zeros(1,M.Nz);
    
    % Advect temperature
    dT(1:end-1) = -dr(1:end-1) .* (M.T(2:end)-M.T(1:end-1))./M.dr;
    M.T  = M.T + dT;
    
    
    %%%%%%%
    % Interpolate new values
    %%%%%%%
    % Volume and mass changes
    dV_s = n2sVolumetric(M,dV_n);
    dm_s = -dV_s*BOD.rhoIce_0;
    dm_s(dm_s<0) = 0;
    
    % Update boundaries
    M.rOcn = interp1(M.r,rTemp,M.rOcn);
    M.rSil = interp1(M.r,rTemp,M.rSil);
    M.rIrn = interp1(M.r,rTemp,M.rIrn);
    M.rIce = interp1(M.r,rTemp,BOD.R);
    
    % Update indices
    M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last'); % Ocean top interface element index
    M.iOcnBot = find((M.rSil - M.r)>=0,1,'last'); % Ocean bottom interface element index
    M.iSilBot = find((M.rIrn - M.r)>=0,1,'last'); % Ocean bottom interface element index
    
    %%%%%%%
    % Handle erupted water
    %%%%%%%
    % Note! For simplicity, it is easy to hangle this as super-heated ice above
    % the melting temperature accreted to tbe base of the ice shell and then
    % fed into the melting/freezing solver to handle appropriately.
    
    if M.iOcnBot==M.iOcnTop
        % Boundaries open up
        M.rOcn = ((3/(4*pi))*(-dV_n(M.iOcnTop)) + M.rOcn^3).^(1/3);
        M.rSil = ((3/(4*pi))*(-dV_n(M.iOcnTop)) + M.rSil^3).^(1/3);
    else % Seafloor will naturally advect down now
        % Boundaries open up
        M.rOcn = ((3/(4*pi))*(-dV_n(M.iOcnTop)) + M.rOcn^3).^(1/3);
    end
    M.vOcn = (4/3)*pi*(M.rOcn.^3 - M.rSil.^3);
    M.vH2O = (4/3)*pi*(BOD.R.^3 - M.rOcn.^3);
    
    % Update indices
    M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last'); % Ocean top interface element index
    M.iOcnBot = find((M.rSil - M.r)>=0,1,'last'); % Ocean bottom interface element index
    M.iSilBot = find((M.rIrn - M.r)>=0,1,'last'); % Ocean bottom interface element index
    
    % Update temperature
    ind = M.iOcnTop;
    T1  = M.T(ind); rho1 = M.rho(ind); Cp1 = M.Cp(ind); V1 = M.V(ind);
    E1  = T1*rho1*Cp1*V1;
    
    T2  = IN.Tm_ocn; rho2 = BOD.rhoIce_0; Cp2 = BOD.CpIce_0; V2 = sum(-dV);
    E2  = T2*rho2*Cp2*V2+sum(dm)*BOD.LH2O;
    
    M.T(ind) = (E1+E2)/((rho1*Cp1*V1+rho2*Cp2*V2));
    
    %%%%%%%
    % Composition info
    %%%%%%%
    % Volume and mass fractions on elements
    fV_H2O_s = zeros(1,M.Nz-1); 
    fV_sil_s = zeros(1,M.Nz-1); 
    fV_irn_s = zeros(1,M.Nz-1); 

    % First assign elements contained within layers
    fV_sil_s(M.iSilBot+1:M.iOcnBot-1) = 1;

    % Now assign boundary elements
    % Seafloor
    fV_sil_s(M.iOcnBot) = (4/3)*pi*(M.rSil^3-M.r(M.iOcnBot)^3)/M.V_s(M.iOcnBot);

    % Core-mantle boundary
    fV_sil_s(M.iSilBot) = (4/3)*pi*(M.r(M.iSilBot+1)^3-M.rIrn^3)/M.V_s(M.iSilBot);

    % Differentiated mass fraction
    dfDiff = (dm_s./(IN.fm0_H2O.*BOD.rhoBulk.*M.V_s))./fV_sil_s;
    dfDiff(isnan(dfDiff) | isinf(dfDiff)) = 1;
    M.fmDiff_H2O = M.fmDiff_H2O+dfDiff;
    M.fmDiff_H2O(M.fmDiff_H2O<0) = 0;
    M.fmDiff_H2O(M.fmDiff_H2O>1) = 1;
    
    % Now interpolate
    rTemp_s = n2sVolumetric(M,rTemp);
    M.fmDiff_H2O(2:end-1) = interp1(rTemp_s,M.fmDiff_H2O,M.r_s(2:end-1));
    
    % Leech K
    M.fmK = IN.fm0_k .* M.fmDiff_H2O;
    
end
end



















