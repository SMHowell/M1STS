%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [M] = getThermalProperties(M,BOD,IN)

%%%%%%%%%%%%%%%%%%
% H2O
%%%%%%%%%%%%%%%%%%
% Conductivity
k0_c   = 2.107; % Ref thermal cond of continuous phase @ BOD.Tice_0 [W/ m K]
k0_d   = 0;     % Ref thermal cond of dispersed phase @ BOD.Tice_0 [W/ m K]
kT_c   = k0_c*BOD.Tice_0./M.T;   % Continuous phase conductivity [W/m K]
kT_d   = k0_d;  % Dispersed phase conductivity [W/m K]
 
k_s_fp = kT_c .* (M.phi - 1).^2 + kT_d; % Snow
k_f_fp = kT_c .* ((1-M.phi) + M.phi .* (3 * k0_d) ./ (2.*kT_c + kT_d))./((1-M.phi) + M.phi .* (3 * kT_c) ./ (2*kT_c + kT_d)); % Firn
k_H2O  = (M.phi./k_s_fp + (1-M.phi)./(k_f_fp)).^-1; % Composite
 
% Heat capacity
Cp_H2O =  1e3 * (7.73e-3 * M.T .* (1 - exp(-1.263e-3 * M.T.^2)) .* ... % Specific heat capacity [J/kg K] 10.1051/0004-6361:20031746
           (1 + exp(-3 * M.T.^(1/2)) * 8.47e-3 .* M.T.^6 + ...
            2.0825e-7 * M.T.^4 .* exp(-4.97e-2 * M.T)));       
    
% Thermal Expansivity [1/K]
a0        = 1.704e-4; % Ref thermal expansivity [1/K]
alpha_ice = a0 * M.T / BOD.Tice_0; % Thermal expansivity [1/K]
 
% Density
rho_H2O   = (1-M.phi).*BOD.rhoIce_0.*(1-(M.T-BOD.Tice_0).*alpha_ice);  % Density [kg/m^3]
        
% Account for melt contribution to Cp, rho, but *NOT* k, which is handled in
% the thermal solver and convectiveConductivity.m. 
% Lay out ocean and reservoir thermal properties.
rhoMelt = BOD.rhoOcn * ones(1,M.Nz);
CpMelt  = BOD.CpOcn * ones(1,M.Nz);
 
% Get Density via volume melt fraction
vfm     = s2nVolumetric(M,M.vfm); % Melt fraction on nodes
rho_H2O = rho_H2O .* (1-vfm) + rhoMelt .* vfm; % Density
 
% Get Cp via mass melt fraction
Nz      = numel(M.r);
rho_s   = n2sVolumetric(M,rho_H2O); % Staggered density
refIndT = 2:Nz-1;                   % Indices of solution
mfm     = zeros(1,Nz);              % Mass melt fraction
mfm(refIndT) = (M.vfm(1:end-1).*M.V_s(1:end-1).*rho_s(1:end-1) + ... 
                M.vfm(2:end).*M.V_s(2:end).*rho_s(2:end))./ ...
               (M.V_s(1:end-1).*rho_s(1:end-1)+M.V_s(2:end).*rho_s(1:end-1));
Cp_H2O  = Cp_H2O .* (1-mfm) + CpMelt .* mfm; % Specific heat capacity
 
 
%%%%%%%%%%%%%%%%%%
% ROCK
%%%%%%%%%%%%%%%%%%
% Conductivity
c0 = 2.47;    % [W m-1 K-1]
c1 = 0.33e-9; % [W m-1 K-1 Pa-1]
c2 = 0.48;    % [N/A]
k_sil = (c0 + c1 * M.P).*(BOD.Tsil_0./M.T).^c2;
 
% Assume forsterite -- Wadsleyite (spinel) transition not until ~15 GPa
% https://www.sciencedirect.com/science/article/pii/S0031920113000289
% Thermal Expansivity [1/K]
a0 = 3.15e-5;  % [K-1]
a1 = 1.02e-8;  % [K-2]
a2 = -0.76;    % [K]
a3 = 3.63e-11; % [Pa-1]
 
alpha_sil = (a0 + a1*M.T + a2*M.T.^(-2)).*exp(-a3*M.P);
 
% Density [kg m3]
rho_sil = (1-M.phi).*BOD.rhoSil_0.*(1-(M.T-BOD.Tsil_0).*alpha_sil); 
 
% Heat capacity 
% High Temps (Above 383 K)
% Gillet, Philippe, et al. "High?temperature thermodynamic properties of 
% forsterite." Journal of Geophysical Research: Solid Earth 96.B7 (1991): 11805-11816
Cp0 = -402.753;
Cp1 = 74.290;
Cp2 = 87.588e3;
Cp3 = -25.913e6;
Cp4 = 25.374e8;
Fo_mm = 140.69/1e3; % Molar mass (kg / mol)
 
Cp_sil_warm = (Cp0 + Cp1*log(M.T) + Cp2./M.T + Cp3./(M.T).^2 + Cp4./(M.T).^3)/Fo_mm;

% Low temps (Below 383K) "Heat capacity and entropy of fayalite (Fe2SiO4)
% between 5.1 and 383 K: comparison of calorimetric and equilibrium values 
% for the QFM buffer reaction" Richard A. Robie; Cabell B. Finch; Bruce S. Hemingway
% American Mineralogist (1982) 67 (5-6): 463–469.
T_cold = 380; % Upper temperature limit of this approx.
Cp0 = 176.02;
Cp1 = -8.808e-3;
Cp2 = 2.471e-5;
Cp3 = -3.889e6;
Fo_mm = 140.69/1e3; % Molar mass (kg / mol)
 
Cp_sil_cold = (Cp0 + Cp1*M.T + Cp2.*M.T.^2 + Cp3.*(M.T).^(-2))/Fo_mm;

% Now for much lower temps, same source
T_veryCold = 180; % Upper temperature limit of this approx.
Cp0  = 4.37e-3;
a_Cp = 2.304;

Cp_sil_veryCold = Cp0*M.T.^a_Cp;

% Hybridize the curves by smoothing with tanh
Sa = T_veryCold; % Smoothing switch point
Sb = T_veryCold/5; % Smoothing window
S_veryCold = 0.5+0.5*tanh((M.T-Sa)/Sb); 

Sa = T_cold; % Smoothing switch point
Sb = T_cold/5; % Smoothing window
S_cold = 0.5+0.5*tanh((M.T-Sa)/Sb); 

Cp_sil = Cp_sil_veryCold .* (1-S_veryCold) + Cp_sil_cold .* (S_veryCold);
Cp_sil = Cp_sil .* (1-S_cold) + Cp_sil_warm .* (S_cold);


%%%%%%%%%%%%%%%%%%
% IRON
%%%%%%%%%%%%%%%%%%
% *********** NOTE: Need to introduce temperature, presure ********
k_irn   = BOD.kIrn_0*ones(1,M.Nz);
rho_irn = BOD.rhoIrn_0*ones(1,M.Nz);
Cp_irn  = BOD.CpIrn_0*ones(1,M.Nz);
 
 
%%%%%%%%%%%%%%%%%%
% Composite
%%%%%%%%%%%%%%%%%%
% Volume and mass fractions on elements
fV_H2O_s = zeros(1,M.Nz-1); 
fV_sil_s = zeros(1,M.Nz-1); 
fV_irn_s = zeros(1,M.Nz-1); 

% First assign elements contained within layers
fV_H2O_s(M.iOcnBot+1:end)         = 1;
fV_sil_s(M.iSilBot+1:M.iOcnBot-1) = 1;
fV_irn_s(1:M.iSilBot-1)           = 1;

% Now assign boundary elements
% Seafloor
fV_H2O_s(M.iOcnBot) = (4/3)*pi*(M.r(M.iOcnBot+1)^3-M.rSil^3)/M.V_s(M.iOcnBot);
fV_sil_s(M.iOcnBot) = (4/3)*pi*(M.rSil^3-M.r(M.iOcnBot)^3)/M.V_s(M.iOcnBot);

% Core-mantle boundary
fV_sil_s(M.iSilBot) = (4/3)*pi*(M.r(M.iSilBot+1)^3-M.rIrn^3)/M.V_s(M.iSilBot);
fV_irn_s(M.iSilBot) = (4/3)*pi*(M.rIrn^3-M.r(M.iSilBot)^3)/M.V_s(M.iSilBot);

% volume fractions on nodes
[fV_H2O] = s2nVolumetric(M,fV_H2O_s);
[fV_sil] = s2nVolumetric(M,fV_sil_s);
[fV_irn] = s2nVolumetric(M,fV_irn_s);

% Mass fractions on nodes
m_H2O = fV_H2O .* rho_H2O .* M.V; % Mass of ice
m_sil = fV_sil .* rho_sil .* M.V; % Mass of rock 
m_irn = fV_irn .* rho_irn .* M.V; % Mass of iron
m_tot = (m_H2O + m_sil + m_irn); % Mass total

fm_H2O = m_H2O ./ m_tot;
fm_sil = m_sil ./ m_tot;
fm_irn = m_irn ./ m_tot;


% Consider differentiation
% Get mass fractions of differentiated material on nodes
[fmDiff_H2O] = s2nVolumetric(M,M.fmDiff_H2O);
[fmDiff_irn] = s2nVolumetric(M,M.fmDiff_irn);

dfm_H2O = fm_sil .* (1-fmDiff_H2O) * IN.fm0_H2O;
dfm_irn = fm_sil .* (1-fmDiff_irn) * IN.fm0_irn;

fm_H2O  = fm_H2O + dfm_H2O;
fm_irn  = fm_irn + dfm_irn;
fm_sil  = fm_sil - (dfm_H2O + dfm_irn);

% Apply to volumes
% Get ice from rock-iron
fm_RI_irn = fm_irn./(fm_irn + fm_sil);
fV_RI_irn = fm_RI_irn./(fm_RI_irn+(1-fm_RI_irn).*rho_irn./rho_sil);
rho_RI    = fV_RI_irn.*rho_irn + (1-fV_RI_irn).*rho_sil;

% Now calculate ice volume fraction!
fV_H2O = fm_H2O./(fm_H2O+(1-fm_H2O).*rho_H2O./rho_RI);
fV_H2O(isnan(fV_H2O)) = 1;

% Rock-Ice portion
fm_RH_H2O = fm_H2O./(fm_H2O + fm_sil);
fV_RH_H2O = fm_RH_H2O./(fm_RH_H2O+(1-fm_RH_H2O).*rho_H2O./rho_sil);
rho_RH    = fV_RH_H2O.*rho_H2O + (1-fV_RH_H2O).*rho_sil;

% Now calculate iron volume fraction!
fV_irn = fm_irn./(fm_irn+(1-fm_irn).*rho_irn./rho_RH);
fV_irn(isnan(fV_irn)) = 1;

% Ice-Iron portion
fm_IH_H2O = fm_H2O./(fm_H2O + fm_irn);
fV_IH_H2O = fm_IH_H2O./(fm_IH_H2O+(1-fm_IH_H2O).*rho_H2O./rho_irn);
rho_IH    = fV_IH_H2O.*rho_H2O + (1-fV_IH_H2O).*rho_irn;

% Now calculate silicate volume fraction!
fV_sil = fm_sil./(fm_sil+(1-fm_sil).*rho_sil./rho_IH);
fV_sil(isnan(fV_sil)) = 1;

% Save fractions *on nodes* 
M.fm_H2O  = fm_H2O;
M.fm_irn  = fm_irn;
M.fm_sil  = fm_sil;

M.fV_H2O  = fV_H2O;
M.fV_irn  = fV_irn;
M.fV_sil  = fV_sil;

% Now mix
M.k   = (fV_H2O./k_H2O  + fV_sil./k_sil   + fV_irn./k_irn).^(-1); 
M.rho = rho_H2O.*fV_H2O + rho_sil.*fV_sil + rho_irn.*fV_irn; 
M.Cp  = Cp_H2O .*fm_H2O + Cp_sil .*fm_sil + Cp_irn .*fm_irn;

% Diffusivity
M.K = M.k./(M.rho .* M.Cp);



end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

