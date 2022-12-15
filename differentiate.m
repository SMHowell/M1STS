%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = differentiate(M,IN,BOD,MAT)

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
fm_n   = s2nVolumetric(M,M.mat.fm_s);
fm_sil = fm_n(M.mat.iSilSolid,:)+fm_n(M.mat.iSilMelt,:);
fm_H2O = fm_n(M.mat.iH2Osolid,:);
ind    = find(fm_H2O>0 & fm_sil>0 & M.T>MAT.H2O.Tm0 & fm_H2O>0);
ind(ind==M.iOcnBot | ind==1) = [];

% Energy available
dE  = 2*pi.*M.rho(ind).*M.Cp(ind).*M.r(ind).^2.*...
    (M.dr(ind-1)+M.dr(ind)).*(M.T(ind)-MAT.H2O.Tm0);

% Total energy required to melt water present
% Get mass fractions of differentiated material on nodes
dE_max = fm_H2O(ind).*M.V(ind).*M.rho(ind)*MAT.H2O.L;

%%%%%%%
% Calculate mass change
%%%%%%%
dExcess = dE-dE_max;
dExcess(dExcess<0) = 0;

% Melted mass of water on node
dm = min(dE,dE_max)./MAT.H2O.L;
if sum(dm)<0; error('Something is broken, you''re trying to bury water'); end

if any(dm>0)
    % Because we withdrew the extra energy, make sure to set
    % temperature accordingly.
    M.T(ind) = MAT.H2O.Tm0;
    
    dT = dExcess./(2*pi*M.rho(ind).*M.Cp(ind).*M.r(ind).^2.*...
        (M.dr(ind-1)+M.dr(ind)));
    M.T(ind) = M.T(ind) + dT;
    
    %%%%%%%
    % Account for mass loss
    %%%%%%%
    % Change in rock volume on node
    dV   = -dm./MAT.H2O.m.rho0;
    dV_n = zeros(1,M.Nz); dV_n(ind) = dV;    
    
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
    M.T = M.T + dT;
    
    %%%%%%%
    % Remove the mass from the rock
    %%%%%%%
    massFull = fm_n.*M.rho.*M.V;
    massFull(M.mat.iH2Osolid,ind) = massFull(M.mat.iH2Osolid,ind) - dm;
    fm_new = massFull./sum(massFull);
    dfm    = fm_new - fm_n;
    dfm_s  = n2sVolumetric(M,dfm);
    M.mat.fm_s = M.mat.fm_s + dfm_s;
    
    rhoFull_s = n2sVolumetric(M,M.mat.rhoFull);
    drho_s  = sum(M.mat.fm_s./rhoFull_s).^-1-M.rho_s;
    drho    = s2nVolumetric(M,drho_s);
    M.rho   = M.rho+drho;
    M.rho_s = M.rho_s+drho_s;
    
    M.mat.fV_s = (M.rho_s .* M.mat.fm_s ./ rhoFull_s); 
    
        
    %%%%%%%
    % Leech Kc (on elements)
    %%%%%%%
    dfK = - IN.fK0 .* dfm_s(M.mat.iH2Osolid,ind);
    dfK(dfK<0) = 0; % Cant give it back!
    M.mat.fK(ind) = M.mat.fK(ind) + dfK;
    M.mat.fK(M.mat.fK>1) = 1; % A limit to all things
    
    
    %%%%%%%
    % Interpolate new values
    %%%%%%%
    % Advect composition 
    dr_s = n2sVolumetric(M,dr);
    
    for i=1:size(M.mat.fm_s,1)
        dfm = interp1(M.r_s+dr_s,M.mat.fm_s(i,:),M.r_s)-M.mat.fm_s(i,:);
        dfV = interp1(M.r_s+dr_s,M.mat.fV_s(i,:),M.r_s)-M.mat.fV_s(i,:);
        
        M.mat.fm_s(i,1:end-1) = M.mat.fm_s(i,1:end-1) + dfm(1:end-1);
        M.mat.fV_s(i,1:end-1) = M.mat.fV_s(i,1:end-1) + dfV(1:end-1);
    end
    dfK = interp1(M.r_s+dr_s,M.mat.fK,M.r_s)-M.mat.fK;
    M.mat.fK(1:end-1) = M.mat.fK(1:end-1) + dfK(1:end-1);
    
    % Make sure nothing got weird during advection
    M.mat.fm_s(M.mat.fm_s<0) = 0;
    M.mat.fm_s(M.mat.fV_s<0) = 0;
    M.mat.fm_s(M.mat.fm_s>1) = 1;
    M.mat.fm_s(M.mat.fV_s>1) = 1;
    
    % Update boundaries
    M.rOcn = interp1(M.r,rTemp,M.rOcn);
    M.rSil = interp1(M.r,rTemp,M.rSil);
    M.rIrn = interp1(M.r,rTemp,M.rIrn);
    
    % Update indices
    M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last'); % Ocean top interface element index
    M.iOcnBot = find((M.rSil - M.r)>=0,1,'last'); % Ocean bottom interface element index
    M.iSilBot = find((M.rIrn - M.r)>=0,1,'last'); % CMB interface element index
    
    %%%%%%%
    % Handle erupted water
    %%%%%%%
    % Note! For simplicity, it is easy to handle this as super-heated ice above
    % the melting temperature accreted to the base of the ice shell and then
    % fed into the melting/freezing solver to handle appropriately.
    M.vOcn = (4/3)*pi*(M.rOcn.^3 - M.rSil.^3);
    M.vH2O = (4/3)*pi*(BOD.R.^3 - M.rOcn.^3);
    
    % Update indices
    M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last'); % Ocean top interface element index
    M.iOcnBot = find((M.rSil - M.r)>=0,1,'last'); % Ocean bottom interface element index
    M.iSilBot = find((M.rIrn - M.r)>=0,1,'last'); % Ocean bottom interface element index
    
    % Update temperature
    ind  = M.iOcnTop;
    fm_n = s2nVolumetric(M,M.mat.fm_s);
    T1   = M.T(ind); rho1 = M.rho(ind); Cp1 = M.Cp(ind); V1 = M.V(ind);
    E1   = T1*rho1*Cp1*V1 + fm_n(M.mat.iH2Omelt,ind).*M.rho(ind).*M.V(ind)*MAT.H2O.L;
    
    T2   = MAT.H2O.Tm0; rho2 = MAT.H2O.s.rho0; Cp2 = MAT.H2O.s.Cp0; V2 = sum(-dV);
    E2   = T2*rho2*Cp2*V2+sum(dm)*MAT.H2O.L;
    
    M.T(ind) = (E1+E2)/((rho1*Cp1*V1+rho2*Cp2*V2));
    
    
end






end



















