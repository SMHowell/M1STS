%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [M] = getThermalProperties(M,IN)
%%%%%%%%%%%%%%%%%%%%%%%
% Initialize thermal properties
%%%%%%%%%%%%%%%%%%%%%%%
rhoOcn   = 1000;  % Density [kg/m^3]
CpOcn    = 4184;  % Specific heat capacity [J/kg K]
LH2O     = 330e3; % Latent heat of fusion [kJ/kg]
rho0_ice = 917;   % Ice reference density [kg/m^3]
T0_ice   = 273;   % Reference temeprature for laws [K]
 
% Rock: https://www.sciencedirect.com/science/article/pii/S0031920113000289
% using forsterite
Tm_sil   = 1300;  % Melting temp of rock [K] 
k_sil    = 4.0;   % Thermal conductivity [W/mK]
Lsil     = 418e3; % Latent heat of fusion [kJ/kg]
rho0_sil = 3275;  % Reference Density [kg/m^3]
T0_sil   = 300;   % Reference temeprature for laws [K] 
 
Tm_irn     = 1260;  % Melting temp of iron [K] 
k0_irn     = 27.0;  % Thermal conductivity [W/mK]
rho0_irn   = 5428;  % Density [kg/m^3]
rho0_irn_m = 5240;  % Density [kg/m^3]
Cp0_irn    = 450;   % Specific heat capacity [J/kg K]
Cp0_irn_m  = 820;   % Specific heat capacity of melt [J/kg K]
Lirn       = 247e3; % Latent heat of fusion [kJ/kg]
 
%%%%%%%%%%%%%%%%%%
% H2O
%%%%%%%%%%%%%%%%%%
% Conductivity
k0_c   = 2.107; % Ref thermal cond of continuous phase @ T0_ice [W/ m K]
k0_d   = 0;     % Ref thermal cond of dispersed phase @ T0_ice [W/ m K]
kT_c   = k0_c*T0_ice./M.T;   % Continuous phase conductivity [W/m K]
kT_d   = k0_d;  % Dispersed phase conductivity [W/m K]
 
k_s_fp = kT_c .* (M.phi - 1).^2 + kT_d; % Snow
k_f_fp = kT_c .* ((1-M.phi) + M.phi .* (3 * k0_d) ./ (2.*kT_c + kT_d))./((1-M.phi) + M.phi .* (3 * kT_c) ./ (2*kT_c + kT_d)); % Firn
k_ice  = (M.phi./k_s_fp + (1-M.phi)./(k_f_fp)).^-1; % Composite
 
% Heat capacity
Cp_ice =  1e3 * (7.73e-3 * M.T .* (1 - exp(-1.263e-3 * M.T.^2)) .* ... % Specific heat capacity [J/kg K] 10.1051/0004-6361:20031746
           (1 + exp(-3 * M.T.^(1/2)) * 8.47e-3 .* M.T.^6 + ...
            2.0825e-7 * M.T.^4 .* exp(-4.97e-2 * M.T)));       
    
% Thermal Expansivity [1/K]
a0        = 1.704e-4; % Ref thermal expansivity [1/K]
alpha_ice = a0 * M.T / T0_ice; % Thermal expansivity [1/K]
 
% Density
rho_ice   = (1-M.phi).*rho0_ice.*(1-(M.T-T0_ice).*alpha_ice);  % Density [kg/m^3]
        
% Account for melt contribution to Cp, rho, but *NOT* k, which is handled in
% the thermal solver and convectiveConductivity.m. 
% Lay out ocean and reservoir thermal properties.
rhoMelt = rhoOcn * ones(1,M.Nz);
CpMelt  = CpOcn * ones(1,M.Nz);
if M.vRes>0; CpMelt(M.r>= M.rResBot & M.r<= M.rResTop) = M.CpRes; end
 
% Get Density via volume melt fraction
vfm     = s2nVolumetric(M,M.vfm); % Melt fraction on nodes
rho_ice = rho_ice .* (1-vfm) + rhoMelt .* vfm; % Density
 
% Get Cp via mass melt fraction
Nz      = numel(M.r);
rho_s   = n2sVolumetric(M,rho_ice); % Staggered density
refIndT = 2:Nz-1;                   % Indices of solution
mfm     = zeros(1,Nz);              % Mass melt fraction
mfm(refIndT) = (M.vfm(1:end-1).*M.V_s(1:end-1).*rho_s(1:end-1) + ... 
                M.vfm(2:end).*M.V_s(2:end).*rho_s(2:end))./ ...
               (M.V_s(1:end-1).*rho_s(1:end-1)+M.V_s(2:end).*rho_s(1:end-1));
Cp_ice  = Cp_ice .* (1-mfm) + CpMelt .* mfm; % Specific heat capacity
 
 
%%%%%%%%%%%%%%%%%%
% ROCK
%%%%%%%%%%%%%%%%%%
% Conductivity
c0 = 2.47;    % [W m-1 K-1]
c1 = 0.33e-9; % [W m-1 K-1 Pa-1]
c2 = 0.48;    % [N/A]
k_sil = (c0 + c1 * M.P).*(T0_sil./M.T).^c2;
 
% Assume forsterite -- Wadsleyite (spinel) transition not until ~15 GPa
% https://www.sciencedirect.com/science/article/pii/S0031920113000289
% Thermal Expansivity [1/K]
a0 = 3.15e-5;  % [K-1]
a1 = 1.02e-8;  % [K-2]
a2 = -0.76;    % [K]
a3 = 3.63e-11; % [Pa-1]
 
alpha_sil = (a0 + a1*M.T + a2*M.T.^(-2)).*exp(-a3*M.P);
 
% Density [kg m3]
rho_sil = (1-M.phi).*rho0_sil.*(1-(M.T-T0_sil).*alpha_sil); 
 
% Heat capacity
% Gillet, Philippe, et al. "High?temperature thermodynamic properties of 
% forsterite." Journal of Geophysical Research: Solid Earth 96.B7 (1991): 11805-11816
Cp0 = -402.753;
Cp1 = 74.290;
Cp2 = 87.588e3;
Cp3 = -25.913e6;
Cp4 = 25.374e8;
Fo_mm = 140.69/1e3; % Molar mass (kg / mol)
 
Cp_sil = (Cp0 + Cp1*log(M.T) + Cp2./M.T + Cp3./(M.T).^2 + Cp4./(M.T).^3)/Fo_mm;
 
% Correct for water present
rho_sil = rho_sil .* (1-vfm) + rhoMelt .* vfm; % Density
rho_s   = n2sVolumetric(M,rho_sil); % Staggered density
refIndT = 2:Nz-1;                   % Indices of solution
mfm     = zeros(1,Nz);              % Mass melt fraction
mfm(refIndT) = (M.vfm(1:end-1).*M.V_s(1:end-1).*rho_s(1:end-1) + ... 
                M.vfm(2:end).*M.V_s(2:end).*rho_s(2:end))./ ...
               (M.V_s(1:end-1).*rho_s(1:end-1)+M.V_s(2:end).*rho_s(1:end-1));
Cp_sil  = Cp_sil  .* (1-mfm) + CpMelt  .* mfm; % Specific heat capacity


%%%%%%%%%%%%%%%%%%
% IRON
%%%%%%%%%%%%%%%%%%
% *********** NOTE: Need to introduce temperature, presure ********
k_irn   = k0_irn*ones(1,M.Nz);
rho_irn = rho0_irn*ones(1,M.Nz);
Cp_irn  = Cp0_irn*ones(1,M.Nz);
 
 
%%%%%%%%%%%%%%%%%%
% Composite
%%%%%%%%%%%%%%%%%%
iceInd = find(M.r>M.rSil); 
silInd = find(M.r>M.rIrn & M.r<=M.rSil); 
irnInd = find(M.r<=M.rIrn); 
 
M.k(iceInd)   = k_ice(iceInd);
M.rho(iceInd) = rho_ice(iceInd);
M.Cp(iceInd)  = Cp_ice(iceInd);
 
M.k(silInd)   = k_sil(silInd);
M.rho(silInd) = rho_sil(silInd);
M.Cp(silInd)  = Cp_sil(silInd);
 
M.k(irnInd)   = k_irn(irnInd);
M.rho(irnInd) = rho_irn(irnInd);
M.Cp(irnInd)  = Cp_irn(irnInd);
 
% Diffusivity
M.K = M.k./(M.rho .* M.Cp);
 
 
%%%%%%%%%%%%%%%%%%
% Stored Values
%%%%%%%%%%%%%%%%%%
M.Tm_ocn = IN.Tm_ocn; % Melt temperature in ocean
M.rhoOcn = rhoOcn;    % Density in ocean
M.CpOcn  = CpOcn;     % Specific heat capacity in ocean
M.LH2O   = LH2O;      % Latent heat of fusion 

M.Tm_res = IN.Tm_ocn; % Melt temperature in reservoir
M.rhoRes = rhoOcn;    % Density in reservoir
M.CpRes  = CpOcn;     % Specific heat capacity in reservoir

M.Lsil   = Lsil;      % Latent heat of fusion 
 
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

