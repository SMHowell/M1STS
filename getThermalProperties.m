%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = getThermalProperties(M,IN,MAT)

%%%%%%%%%%%%%%%%%%
% H2O Solid
%%%%%%%%%%%%%%%%%%
% Conductivity
k0_c   = 2.107; % Ref thermal cond of continuous phase @ MAT.H2O.s.T0 [W/ m K]
k0_d   = 0;     % Ref thermal cond of dispersed phase @ MAT.H2O.s.T0 [W/ m K]
kT_c   = k0_c*MAT.H2O.s.T0./M.T;   % Continuous phase conductivity [W/m K]
kT_d   = k0_d;  % Dispersed phase conductivity [W/m K]

k_s_fp = kT_c .* (M.phi - 1).^2 + kT_d; % Snow
k_f_fp = kT_c .* ((1-M.phi) + M.phi .* (3 * k0_d) ./ (2.*kT_c + kT_d))./((1-M.phi) + M.phi .* (3 * kT_c) ./ (2*kT_c + kT_d)); % Firn
kH2O_s = (M.phi./k_s_fp + (1-M.phi)./(k_f_fp)).^-1; % Composite

% Heat capacity
CpH2O_s =  1e3 * (7.73e-3 * M.T .* (1 - exp(-1.263e-3 * M.T.^2)) .* ... % Specific heat capacity [J/kg K] 10.1051/0004-6361:20031746
    (1 + exp(-3 * M.T.^(1/2)) * 8.47e-3 .* M.T.^6 + ...
    2.0825e-7 * M.T.^4 .* exp(-4.97e-2 * M.T)));

% Thermal Expansivity [1/K]
a0        = 1.704e-4; % Ref thermal expansivity [1/K]
alpha_ice = a0 * M.T / MAT.H2O.s.T0; % Thermal expansivity [1/K]

% Density
rhoH2O_s = (1-M.phi).*MAT.H2O.s.rho0.*(1-(M.T-MAT.H2O.s.T0).*alpha_ice);  % Density [kg/m^3]


%%%%%%%%%%%%%%%%%%
% H2O Melt
%%%%%%%%%%%%%%%%%%
% Account for melt contribution to Cp, rho, but *NOT* k, which is handled in
% the thermal solver and convectiveConductivity.m.
% Lay out ocean and reservoir thermal properties.
rhoH2O_m = MAT.H2O.m.rho0 * ones(1,M.Nz);
CpH2O_m  = MAT.H2O.m.Cp0  * ones(1,M.Nz);


%%%%%%%%%%%%%%%%%%
% Rock
%%%%%%%%%%%%%%%%%%
switch IN.silComp
    
    case 'fo'
        %%%%%%%%%%%%%%%%%%
        % Forsterite
        %%%%%%%%%%%%%%%%%%
        % Melting Tempertature
        % Get volume change fV = V/V0
        K  = MAT.SIL.s.Kmod; % Bulk modulus
        fv = exp(-M.P./K);
        C  = 3;
        M.TmSil = MAT.SIL.Tm0 * (1 + C * (1-fv));
        
        % Conductivity
        c0 = 2.47;    % [W m-1 K-1]
        c1 = 0.33e-9; % [W m-1 K-1 Pa-1]
        c2 = 0.48;    % [N/A]
        kSil_s = (c0 + c1 * M.P).*(MAT.SIL.s.T0./M.T).^c2;
        
        % Assume forsterite -- Wadsleyite (spinel) transition not until ~15 GPa
        % https://www.sciencedirect.com/science/article/pii/S0031920113000289
        % Thermal Expansivity [1/K]
        a0 = 3.15e-5;  % [K-1]
        a1 = 1.02e-8;  % [K-2]
        a2 = -0.76;    % [K]
        a3 = 3.63e-11; % [Pa-1]
        
        alpha_sil = (a0 + a1*M.T + a2*M.T.^(-2)).*exp(-a3*M.P);
        alpha_sil(alpha_sil<1e-6) = 1e-6; % Extrap blows up when cold
        
        % Density [kg m3]
        rhoSil_s = (1-M.phi).*MAT.SIL.s.rho0.*(1-(M.T-MAT.SIL.s.T0).*alpha_sil);
        
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
        
        CpSil_warm = (Cp0 + Cp1*log(M.T) + Cp2./M.T + Cp3./(M.T).^2 + Cp4./(M.T).^3)/Fo_mm;
        
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
        
        CpSil_cold = (Cp0 + Cp1*M.T + Cp2.*M.T.^2 + Cp3.*(M.T).^(-2))/Fo_mm;
        
        % Now for much lower temps, same source
        T_veryCold = 180; % Upper temperature limit of this approx.
        Cp0  = 4.37e-3;
        a_Cp = 2.304;
        
        CpSil_veryCold = Cp0*M.T.^a_Cp;
        
        % Hybridize the curves by smoothing with tanh
        Sa = T_veryCold; % Smoothing switch point
        Sb = T_veryCold/5; % Smoothing window
        S_veryCold = 0.5+0.5*tanh((M.T-Sa)/Sb);
        
        Sa = T_cold; % Smoothing switch point
        Sb = T_cold/5; % Smoothing window
        S_cold = 0.5+0.5*tanh((M.T-Sa)/Sb);
        
        CpSil_s = CpSil_veryCold .* (1-S_veryCold) + CpSil_cold .* (S_veryCold);
        CpSil_s = CpSil_s .* (1-S_cold) + CpSil_warm .* (S_cold);
        
        
    case 'fa'
        %%%%%%%%%%%%%%%%%%
        % Fayalite
        %%%%%%%%%%%%%%%%%%
        % Melting Tempertature
        M.TmSil = MAT.SIL.Tm0 + 4.85e-8 * M.P;
        
        % Conductivity
        c0 = 2.47;    % [W m-1 K-1]
        c1 = 0.33e-9; % [W m-1 K-1 Pa-1]
        c2 = 0.48;    % [N/A]
        kSil_s = (c0 + c1 * M.P).*(MAT.SIL.s.T0./M.T).^c2;
        
        % Assume forsterite -- Wadsleyite (spinel) transition not until ~15 GPa
        % https://www.sciencedirect.com/science/article/pii/S0031920113000289
        % Thermal Expansivity [1/K]
        a0 = 3.15e-5;  % [K-1]
        a1 = 1.02e-8;  % [K-2]
        a2 = -0.76;    % [K]
        a3 = 3.63e-11; % [Pa-1]
        
        alpha_sil = (a0 + a1*M.T + a2*M.T.^(-2)).*exp(-a3*M.P);
        alpha_sil(alpha_sil<1e-6) = 1e-6; % Extrap blows up when cold
        
        % Density [kg m3]
        rhoSil_s = (1-M.phi).*MAT.SIL.s.rho0.*(1-(M.T-MAT.SIL.s.T0).*alpha_sil);
        
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
        
        CpSil_warm = (Cp0 + Cp1*log(M.T) + Cp2./M.T + Cp3./(M.T).^2 + Cp4./(M.T).^3)/Fo_mm;
        
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
        
        CpSil_cold = (Cp0 + Cp1*M.T + Cp2.*M.T.^2 + Cp3.*(M.T).^(-2))/Fo_mm;
        
        % Now for much lower temps, same source
        T_veryCold = 180; % Upper temperature limit of this approx.
        Cp0  = 4.37e-3;
        a_Cp = 2.304;
        
        CpSil_veryCold = Cp0*M.T.^a_Cp;
        
        % Hybridize the curves by smoothing with tanh
        Sa = T_veryCold; % Smoothing switch point
        Sb = T_veryCold/5; % Smoothing window
        S_veryCold = 0.5+0.5*tanh((M.T-Sa)/Sb);
        
        Sa = T_cold; % Smoothing switch point
        Sb = T_cold/5; % Smoothing window
        S_cold = 0.5+0.5*tanh((M.T-Sa)/Sb);
        
        CpSil_s = CpSil_veryCold .* (1-S_veryCold) + CpSil_cold .* (S_veryCold);
        CpSil_s = CpSil_s .* (1-S_cold) + CpSil_warm .* (S_cold);
        
    case 'MB'
        %%%%%%%%%%%%%%%%%%
        % MORB
        %%%%%%%%%%%%%%%%%%
        % Melting Tempertature
        M.TmSil = MAT.SIL.Tm0 + 4.85e-8 * M.P;
        
        % Conductivity
        % https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/JB075i020p04063
        % k = A + B * T^3
        c0 = MAT.SIL.s.k0; % Noninally 2.47 [W m-1 K-1]
        c1 = 0.33e-9; % [W m-1 K-1 Pa-1]
        c2 = 0.48;    % [N/A]
        kSil_s = (c0 + c1 * M.P).*(MAT.SIL.s.T0./M.T).^c2;
        
        % Assume forsterite -- Wadsleyite (spinel) transition not until ~15 GPa
        % https://www.sciencedirect.com/science/article/pii/S0031920113000289
        % Thermal Expansivity [1/K]
        a0 = 3.15e-5;  % [K-1]
        a1 = 1.02e-8;  % [K-2]
        a2 = -0.76;    % [K]
        a3 = 3.63e-11; % [Pa-1]
        
        alpha_sil = (a0 + a1*M.T + a2*M.T.^(-2)).*exp(-a3*M.P);
        alpha_sil(alpha_sil<1e-6) = 1e-6; % Extrap blows up when cold
        
        % Density [kg m3]
        rhoSil_s = (1-M.phi).*MAT.SIL.s.rho0.*(1-(M.T-MAT.SIL.s.T0).*alpha_sil);
        
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
        
        CpSil_warm = (Cp0 + Cp1*log(M.T) + Cp2./M.T + Cp3./(M.T).^2 + Cp4./(M.T).^3)/Fo_mm;
        
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
        
        CpSil_cold = (Cp0 + Cp1*M.T + Cp2.*M.T.^2 + Cp3.*(M.T).^(-2))/Fo_mm;
        
        % Now for much lower temps, same source
        T_veryCold = 180; % Upper temperature limit of this approx.
        Cp0  = 4.37e-3;
        a_Cp = 2.304;
        
        CpSil_veryCold = Cp0*M.T.^a_Cp;
        
        % Hybridize the curves by smoothing with tanh
        Sa = T_veryCold; % Smoothing switch point
        Sb = T_veryCold/5; % Smoothing window
        S_veryCold = 0.5+0.5*tanh((M.T-Sa)/Sb);
        
        Sa = T_cold; % Smoothing switch point
        Sb = T_cold/5; % Smoothing window
        S_cold = 0.5+0.5*tanh((M.T-Sa)/Sb);
        
        CpSil_s = CpSil_veryCold .* (1-S_veryCold) + CpSil_cold .* (S_veryCold);
        CpSil_s = CpSil_s .* (1-S_cold) + CpSil_warm .* (S_cold);
end


%%%%%%%%%%%%%%%%%%
% Silicate Melt
%%%%%%%%%%%%%%%%%%
% Account for melt contribution to Cp, rho, but *NOT* k, which is handled in
% the thermal solver and convectiveConductivity.m.
% Lay out ocean and reservoir thermal properties.
rhoSil_m = MAT.SIL.m.rho0 * ones(1,M.Nz);
CpSil_m  = MAT.SIL.m.Cp0  * ones(1,M.Nz);


%%%%%%%%%%%%%%%%%%
% IRON
%%%%%%%%%%%%%%%%%%
% *********** NOTE: Need to introduce temperature, presure ********
kIrn_s   = MAT.IRN.s.k0*ones(1,M.Nz);
rhoIrn_s = MAT.IRN.s.rho0*ones(1,M.Nz);
CpIrn_s  = MAT.IRN.s.Cp0*ones(1,M.Nz);


%%%%%%%%%%%%%%%%%%
% Iron Melt
%%%%%%%%%%%%%%%%%%
% Account for melt contribution to Cp, rho, but *NOT* k, which is handled in
% the thermal solver and convectiveConductivity.m.
% Lay out ocean and reservoir thermal properties.
rhoIrn_m = MAT.IRN.m.rho0 * ones(1,M.Nz);
CpIrn_m  = MAT.IRN.m.Cp0  * ones(1,M.Nz);


%%%%%%%%%%%%%%%%%%
% Composite
%%%%%%%%%%%%%%%%%%

% Save arryas
M.mat.rhoFull(M.mat.iH2Osolid,:) = rhoH2O_s;
M.mat.rhoFull(M.mat.iH2Omelts,:) = rhoH2O_m;
M.mat.rhoFull(M.mat.iSilSolid,:) = rhoSil_s;
M.mat.rhoFull(M.mat.iSilMelts,:) = rhoSil_m;
M.mat.rhoFull(M.mat.iIrnSolid,:) = rhoIrn_s;
M.mat.rhoFull(M.mat.iIrnMelts,:) = rhoIrn_m;

kFull(M.mat.iH2Osolid,:) = kH2O_s;
kFull(M.mat.iH2Omelts,:) = kH2O_s;
kFull(M.mat.iSilSolid,:) = kSil_s;
kFull(M.mat.iSilMelts,:) = kSil_s;
kFull(M.mat.iIrnSolid,:) = kIrn_s;
kFull(M.mat.iIrnMelts,:) = kIrn_s;

CpFull(M.mat.iH2Osolid,:) = CpH2O_s;
CpFull(M.mat.iH2Omelts,:) = CpH2O_m;
CpFull(M.mat.iSilSolid,:) = CpSil_s;
CpFull(M.mat.iSilMelts,:) = CpSil_m;
CpFull(M.mat.iIrnSolid,:) = CpIrn_s;
CpFull(M.mat.iIrnMelts,:) = CpIrn_m;


% Update density
rhoFull_s  = n2sVolumetric(M,M.mat.rhoFull);
M.rho_s    = sum(M.mat.fm_s./rhoFull_s).^-1;
M.rho      = s2nVolumetric(M,M.rho_s);

% volume fractions and interpolate to nodes
M.mat.fV_s = M.rho_s .* M.mat.fm_s ./ rhoFull_s; 
fV_n       = s2nVolumetric(M,M.mat.fV_s);
fm_n       = s2nVolumetric(M,M.mat.fV_s);

% Now mix
M.k   = sum(fV_n./kFull).^(-1);
M.Cp  = sum(CpFull.*fm_n);

% Diffusivity
M.K = M.k./(M.rho .* M.Cp);



end










































