%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M,OUT] = differentiate(M,IN,BOD,MAT,OUT)

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
%%%%%%%
fm_n   = s2nVolumetric(M,M.mat.fm_s);
fV_n   = s2nMass(M,M.mat.fV_s);
fm_sil = fm_n(M.mat.iSilSolid,:)+fm_n(M.mat.iSilMelt,:);
fm_H2O = fm_n(M.mat.iH2Osolid,:);
fm_ocn = fm_n(M.mat.iH2Omelt,:);
ind    = find(fm_H2O>0 & fm_ocn==0 & fm_sil>0 & M.r<M.rSil & M.T>MAT.H2O.Tm0);
ind(ind==1) = [];

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
removeInd = find(dm<=0);
ind(removeInd)     = [];
dm(removeInd)      = [];
dExcess(removeInd) = [];

if any(dm>0) && M.diffH2O
    % Track progress
    if isempty(OUT.diffH2Ostart); OUT.diffH2Ostart = M.t; 
    end
    
    % Because we withdrew the extra energy, make sure to set
    % temperature accordingly.
    M.T(ind) = MAT.H2O.Tm0;
    
    dT = dExcess./(2*pi*M.rho(ind).*M.Cp(ind).*M.r(ind).^2.*...
        (M.dr(ind-1)+M.dr(ind)));
    M.T(ind) = M.T(ind) + dT;
    
    %%%%%%%
    % Account for mass loss
    %%%%%%%
    % Change in volume on node
    dV   = -dm./M.mat.rhoFull(M.mat.iH2Osolid,ind);
    dV_n = zeros(1,M.Nz); dV_n(ind) = dV;
    
    %%%%%%%
    % Fix fractions
    %%%%%%%
    V_n = fV_n .* M.V;
    V_n(M.mat.iH2Osolid,ind) = V_n(M.mat.iH2Osolid,ind) + dV;
    
    m_n = fm_n .* M.V .* M.rho;
    m_n(M.mat.iH2Osolid,ind) = m_n(M.mat.iH2Osolid,ind) - dm;
    
    fV_new = V_n ./ sum(V_n);
    fm_new = m_n ./ sum(m_n);
    
    dfm_n  = fm_new - fm_n;
    dfV_n  = fV_new - fV_n;
    
    dfm_pass = zeros(size(dfm_n));
    dfm_pass(:,ind) = dfm_n(:,ind);
    
    dfV_pass = zeros(size(dfV_n));
    dfV_pass(:,ind) = dfV_n(:,ind);
    
    dfm_s  = n2sVolumetric(M,dfm_pass);
    dfV_s  = n2sVolumetric(M,dfV_pass);
    
    M.mat.fV_s = M.mat.fV_s + dfV_s;
    M.mat.fV_s(M.mat.fV_s<0) = 0;
    M.mat.fV_s = M.mat.fV_s ./ sum(M.mat.fV_s);
    
    M.mat.fm_s = M.mat.fm_s + dfm_s;
    M.mat.fm_s(M.mat.fm_s<0) = 0;
    M.mat.fm_s = M.mat.fm_s ./ sum(M.mat.fm_s);
    
        
    % Propogate change
    rTemp  = M.r;
    for i=2:M.Nz
        rTemp(i) = ((3/(4*pi))*(dV_n(i)+M.V_s(i-1)) + rTemp(i-1)^3).^(1/3);
    end
    
    % Temperature update using upwind finite differences
    dr = rTemp-M.r;
    dT = zeros(1,M.Nz);
    
    % Advect temperature
    dT(1:end-1) = -dr(1:end-1) .* (M.T(2:end)-M.T(1:end-1))./M.dr;
    M.T = M.T + dT;

    %%%%%%%
    % Interpolate new values
    %%%%%%%
    % Advect composition
    dr_s = n2sVolumetric(M,dr);
    
    for i=1:size(M.mat.fm_s,1)
        dfm = interp1(M.r_s+dr_s,M.mat.fm_s(i,:),M.r_s)-M.mat.fm_s(i,:);
        dfV = interp1(M.r_s+dr_s,M.mat.fV_s(i,:),M.r_s)-M.mat.fV_s(i,:);
        
        dfm(any(isnan(dfm))) = 0;
        dfV(any(isnan(dfV))) = 0;
        
        M.mat.fm_s(i,1:end-1) = M.mat.fm_s(i,1:end-1) + dfm(1:end-1);
        M.mat.fV_s(i,1:end-1) = M.mat.fV_s(i,1:end-1) + dfV(1:end-1);
    end
    
    %%%%%%%
    % Check when to stop and leech Kc 
    %%%%%%%
    m0 =  (IN.fm0_H2O .* BOD.m);
    M.mat.dmH2O_diff = M.mat.dmH2O_diff + sum(dm);
    M.mat.fmH2O_diff = M.mat.dmH2O_diff ./ (IN.fm0_H2O .* BOD.m);
    
    if M.mat.dmH2O_diff>m0
        M.mat.dmH2O_diff = m0;
        M.mat.fmH2O_diff = 1;
        M.mat.fm_s(M.mat.iH2Osolid,M.r_s<M.rOcn) = 0;
        M.mat.fV_s(M.mat.iH2Osolid,M.r_s<M.rOcn) = 0;
        M.diffH2O = 0;
        OUT.diffH2Ostop = M.t;
    end

    % Update boundaries
    M.rOcn = interp1(M.r,rTemp,M.rOcn);
    M.rSil = interp1(M.r,rTemp,M.rSil);
    M.rIrn = interp1(M.r,rTemp,M.rIrn);
    M.rIoc = interp1(M.r,rTemp,M.rIoc);
    
    % Update indices
    M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last'); % Ocean top interface element index
    M.iOcnBot = find((M.rSil - M.r)>=0,1,'last'); % Ocean bottom interface element index
    M.iIocTop = find((M.rIrn - M.r)>=0,1,'last'); % CMB interface element index
    M.iIocBot = find((M.rIoc - M.r)>=0,1,'last'); % OC top interface element index
    
        
    % Make sure nothing got weird during advection
    fm_thresh = 1e-4; % Threshold fraction where differentiation is considered complete
    M.mat.fm_s(M.mat.iSilSolid,M.r_s>M.rSil) = 0;
    M.mat.fm_s(M.mat.iSilMelt,M.r_s>M.rSil)  = 0;
    M.mat.fm_s(M.mat.iIrnSolid,M.r_s>M.rSil) = 0;
    M.mat.fm_s(M.mat.iIrnMelt,M.r_s>M.rSil)  = 0;
    
 
    M.mat.fm_s(M.mat.fm_s<fm_thresh) = 0;
    M.mat.fV_s(M.mat.fV_s<fm_thresh) = 0;
    M.mat.fm_s(M.mat.fm_s>1)    = 1;
    M.mat.fV_s(M.mat.fV_s>1)    = 1;
    M.mat.fm_s = M.mat.fm_s ./ sum(M.mat.fm_s);
    M.mat.fV_s = M.mat.fV_s ./ sum(M.mat.fV_s);
    
    %%%%%%%
    % Handle erupted water
    %%%%%%%
    % Note! For simplicity, it is easy to handle this as super-heated ice above
    % the melting temperature accreted to the base of the ice shell and then
    % fed into the melting/freezing solver to handle appropriately.
    M.vOcn = (4/3)*pi*(M.rOcn.^3 - M.rSil.^3);
    M.vH2O = (4/3)*pi*(BOD.R.^3 - M.rOcn.^3);
    M.vIrn = (4/3)*pi*(M.rIrn^3);          % Core volume
    M.vIoc = (4/3)*pi*(M.rIrn^3-M.rIoc^3); % Molten outer core volume
    
    % Update indices
    M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last'); % Ocean top interface element index
    M.iOcnBot = find((M.rSil - M.r)>=0,1,'last'); % Ocean bottom interface element index
    M.iIocTop = find((M.rIrn - M.r)>=0,1,'last'); % Ocean bottom interface element index
    M.iIocBot = find((M.rIoc - M.r)>=0,1,'last'); % OC top interface element index
    
    % Update temperature
    ind  = M.iOcnBot;
    fm_n = s2nVolumetric(M,M.mat.fm_s);
    T1   = M.T(ind); rho1 = M.rho(ind); Cp1 = M.Cp(ind); V1 = M.V(ind);
    E1   = T1*rho1*Cp1*V1 + fm_n(M.mat.iH2Omelt,ind).*M.rho(ind).*M.V(ind)*MAT.H2O.L;
    
    T2   = MAT.H2O.Tm0; rho2 = MAT.H2O.s.rho0; Cp2 = MAT.H2O.s.Cp0; V2 = sum(-dV);
    E2   = T2*rho2*Cp2*V2+sum(dm)*MAT.H2O.L;
    
    M.T(ind) = (E1+E2)/((rho1*Cp1*V1+rho2*Cp2*V2));
    
    
end





% %%%%%%%%%%%%%%%%%%
% % Extract Iron
% %%%%%%%%%%%%%%%%%%
% % Same thing, but down
% 
% %%%%%%%
% % Set up problem
% %%%%%%%
% % Indices where icy differentiation should occur
% fm_n   = s2nVolumetric(M,M.mat.fm_s);
% fm_sil = fm_n(M.mat.iSilSolid,:)+fm_n(M.mat.iSilMelt,:);
% fm_irn = fm_n(M.mat.iIrnSolid,:);
% ind    = find(fm_irn>0 & fm_sil>0 & M.T>MAT.IRN.Tm0 & fm_irn>0);
% ind(ind==1) = [];
% 
% % Energy available
% dE  = 2*pi.*M.rho(ind).*M.Cp(ind).*M.r(ind).^2.*...
%     (M.dr(ind-1)+M.dr(ind)).*(M.T(ind)-MAT.IRN.Tm0);
% 
% % Total energy required to melt water present
% % Get mass fractions of differentiated material on nodes
% dE_max = fm_irn(ind).*M.V(ind).*M.rho(ind)*MAT.IRN.L;
% 
% %%%%%%%
% % Calculate mass change
% %%%%%%%
% dExcess = dE-dE_max;
% dExcess(dExcess<0) = 0;
% 
% % Melted mass of water on node
% dm = min(dE,dE_max)./MAT.IRN.L;
% if sum(dm)<0; error('Something is broken, you''re trying to bury water');
% end
% 
% if any(dm>0)
%     % Track progress
%     if isempty(OUT.diffIrnStart); OUT.diffIrnStart = M.t; 
%     end
%     
%     % Because we withdrew the extra energy, make sure to set
%     % temperature accordingly.
%     M.T(ind) = MAT.IRN.Tm0;
%     
%     dT = dExcess./(2*pi*M.rho(ind).*M.Cp(ind).*M.r(ind).^2.*...
%         (M.dr(ind-1)+M.dr(ind)));
%     M.T(ind) = M.T(ind) + dT;
%     
%     %%%%%%%
%     % Account for mass loss
%     %%%%%%%
%     % Change in rock volume on node
%     dm   = [-sum(dm),dm];
%     ind  = [M.iIocTop,ind];
%     dV   = -dm./MAT.IRN.m.rho0;
%     dV_n = zeros(1,M.Nz); dV_n(ind) = dV;
%     
%     % Propagate change
%     rTemp  = M.r;
%     if M.iIocTop==1
%         rTemp(1) = ((3/(4*pi))*dV_n(1)).^(1/3);
%     end
%     for i=2:M.Nz
%         rTemp(i) = ((3/(4*pi))*(dV_n(i)+M.V_s(i-1)) + rTemp(i-1)^3).^(1/3);
%     end
%     
%     % Temperature update using upwind finite differences
%     dr = rTemp-M.r;
%     dT = zeros(1,M.Nz);
%     
%     % Advect temperature
%     dT(1:end-1) = -dr(1:end-1) .* (M.T(2:end)-M.T(1:end-1))./M.dr;
%     M.T = M.T + dT;
%     
%     %%%%%%%
%     % Remove the mass from the rock
%     %%%%%%%
%     massFull = fm_n.*M.rho.*M.V;
%     massFull(M.mat.iIrnSolid,ind) = massFull(M.mat.iIrnSolid,ind) - dm;
%     fm_new = massFull./sum(massFull);
%     dfm    = fm_new - fm_n;
%     dfm(any(isnan(dfm))) = [];
%     dfm_s  = n2sVolumetric(M,dfm);
%     M.mat.fm_s = M.mat.fm_s + dfm_s;
%     
%     rhoFull_s = n2sVolumetric(M,M.mat.rhoFull);
%     drho_s  = sum(M.mat.fm_s./rhoFull_s).^-1-M.rho_s;
%     drho    = s2nVolumetric(M,drho_s);
%     M.rho   = M.rho+drho;
%     M.rho_s = M.rho_s+drho_s;
%     
%     M.mat.fV_s = (M.rho_s .* M.mat.fm_s ./ rhoFull_s);
%     
%     
%     %%%%%%%
%     % Interpolate new values
%     %%%%%%%
%     % Advect composition
%     dr_s = n2sVolumetric(M,dr);
%     
%     for i=1:size(M.mat.fm_s,1)
%         dfm = interp1(M.r_s+dr_s,M.mat.fm_s(i,:),M.r_s)-M.mat.fm_s(i,:);
%         dfV = interp1(M.r_s+dr_s,M.mat.fV_s(i,:),M.r_s)-M.mat.fV_s(i,:);
%         
%         dfm(any(isnan(dfm))) = 0;
%         dfV(any(isnan(dfV))) = 0;
%         
%         M.mat.fm_s(i,1:end-1) = M.mat.fm_s(i,1:end-1) + dfm(1:end-1);
%         M.mat.fV_s(i,1:end-1) = M.mat.fV_s(i,1:end-1) + dfV(1:end-1);
%     end
%     dfK = interp1(M.r_s+dr_s,M.mat.fK,M.r_s)-M.mat.fK;
%     dfK(any(isnan(dfm))) = 0;
%     M.mat.fK(1:end-1) = M.mat.fK(1:end-1) + dfK(1:end-1);
%     
%     % Make sure nothing got weird during advection
%     M.mat.fm_s(M.mat.fm_s<0) = 0;
%     M.mat.fm_s(M.mat.fV_s<0) = 0;
%     M.mat.fm_s(M.mat.fm_s>1) = 1;
%     M.mat.fm_s(M.mat.fV_s>1) = 1;
%     
%     % Update boundaries
%     M.rOcn = interp1(M.r,rTemp,M.rOcn);
%     M.rSil = interp1(M.r,rTemp,M.rSil);
%     M.rIrn = interp1(M.r,rTemp,M.rIrn);
%     M.rIoc = interp1(M.r,rTemp,M.rIoc);
%     
%     % Update indices
%     M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last'); % Ocean top interface element index
%     M.iOcnBot = find((M.rSil - M.r)>=0,1,'last'); % Ocean bottom interface element index
%     M.iIocTop = find((M.rIrn - M.r)>=0,1,'last'); % CMB interface element index
%     M.iIocBot = find((M.rIoc - M.r)>=0,1,'last'); % OC top interface element index
%     
%     %%%%%%%
%     % Handle coalesced iron
%     %%%%%%%
%     % Note! For simplicity, it is easy to handle this as super-heated iron above
%     % the melting temperature accreted to the core, and let
%     % melting/freezing handle it.
%     M.vOcn = (4/3)*pi*(M.rOcn.^3 - M.rSil.^3);
%     M.vH2O = (4/3)*pi*(BOD.R.^3 - M.rOcn.^3);
%     M.vIrn = (4/3)*pi*(M.rIrn^3);         % Core volume
%     M.vIoc = (4/3)*pi*(M.rIrn^3-M.rIoc^3);   % Molten outer core volume
%     
%     % Update indices
%     M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last'); % Ocean top interface element index
%     M.iOcnBot = find((M.rSil - M.r)>=0,1,'last'); % Ocean bottom interface element index
%     M.iIocTop = find((M.rIrn - M.r)>=0,1,'last'); % Iron top interface element index
%     M.iIocBot = find((M.rIoc - M.r)>=0,1,'last'); % OC top interface element index
%     
%     % Update temperature
%     ind  = M.iIocTop;
%     fm_n = s2nVolumetric(M,M.mat.fm_s);
%     T1   = M.T(ind); rho1 = M.rho(ind); Cp1 = M.Cp(ind); V1 = M.V(ind);
%     E1   = T1*rho1*Cp1*V1 + fm_n(M.mat.iIrnMelt,ind).*M.rho(ind).*M.V(ind)*MAT.IRN.L;
%     
%     T2   = MAT.IRN.Tm0; rho2 = MAT.IRN.s.rho0; Cp2 = MAT.IRN.s.Cp0; V2 = dV(1);
%     E2   = T2*rho2*Cp2*V2+sum(dm)*MAT.IRN.L;
%     
%     M.T(ind) = (E1+E2)/((rho1*Cp1*V1+rho2*Cp2*V2));
%     
%     
% elseif ~isempty(OUT.diffIrnStart) && ~any(M.T(M.r<M.rSil)<MAT.IRN.Tm0)
%     OUT.diffIrnStop = M.t; 
% end


end



















