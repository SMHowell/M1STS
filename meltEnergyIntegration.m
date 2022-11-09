%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dE = meltEnergyIntegration(r1,r2,T1,T2,T1_corr,T2_corr,rOcnTop,N)
% NOTE! All of this should be done assuming solid ice thermal properties,
% even for water. The thermal solver doesnt care about rho, Cp in melt,
% only k, and that is handled already. The appropraite melt values are used
% for melting/freezing elsewhere.

%% Calculate excess energy in partially melted layer
r   = linspace(r1,r2,N); % Temporary radial array
r_s = sqrt(r(1:end-1) .* r(2:end));    % Staggered locations
V_s = (4/3) * pi * (r(2:end).^3-r(1:end-1).^3); % Volumes
T_err = T2 * (T1/T2).^((r2-r_s)/(r2-r1));  % Erroneous Temperatures
T_cor = T2_corr * (T1_corr/T2_corr).^((r2-r_s)/(r2-rOcnTop)); % Correct Temperatures
T_cor(T_cor>T1_corr) = T1_corr;

% Thermal props
T0    = 270;      % Ref temp for laws [K]
a0    = 1.704e-4; % Ref thermal expansivity [1/K]
rho0  = 917;      % Ref density @ T0 [kg/m^3]

% Erroneous values
alpha_err = a0 * T_err / T0;                 % Thermal expansivity [1/K]
rho_err   = rho0.*(1-(T_err-T0).*alpha_err);  % Density [kg/m^3]
Cp_err    =  1e3 * (7.73e-3 * T_err .* (1 - exp(-1.263e-3 * T_err.^2)) .* ... % Specific heat capacity [J/kg K] 10.1051/0004-6361:20031746
            (1 + exp(-3 * T_err.^(1/2)) * 8.47e-3 .* T_err.^6 + ...
            2.0825e-7 * T_err.^4 .* exp(-4.97e-2 * T_err)));

% Correct Values
alpha_cor = a0 * T_cor / T0;                 % Thermal expansivity [1/K]
rho_cor   = rho0.*(1-(T_cor-T0).*alpha_cor);  % Density [kg/m^3]
Cp_cor    =  1e3 * (7.73e-3 * T_cor .* (1 - exp(-1.263e-3 * T_cor.^2)) .* ... % Specific heat capacity [J/kg K] 10.1051/0004-6361:20031746
            (1 + exp(-3 * T_cor.^(1/2)) * 8.47e-3 .* T_cor.^6 + ...
            2.0825e-7 * T_cor.^4 .* exp(-4.97e-2 * T_cor)));        
        
% Integrate to get excess energy
Eice_err = sum(rho_err.*Cp_err.*T_err.*V_s);
Eice_cor = sum(rho_cor.*Cp_cor.*T_cor.*V_s);
dE       = Eice_err-Eice_cor;


end



















