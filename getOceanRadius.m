%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,rOcn_new] = getOceanRadius(M,BOD,IN,dE)

%%%%%%%%%%%%%%%%%%
% Ice Thermal Properties
%%%%%%%%%%%%%%%%%%
% Get props at appropraite temps
% Handle degenerate case where interface is already at melting temp
if M.T(M.iOcnTop+1) < IN.Tm_ocn
    T_avg = (IN.Tm_ocn-M.T(M.iOcnTop+1))/log(IN.Tm_ocn/M.T(M.iOcnTop+1));
else
    T_avg = IN.Tm_ocn;
end
T   = [IN.Tm_ocn, T_avg];
phi = 0; % Zero out porosity dependence for interface problem

% Heat capacity
Cp_ice =  1e3 * (7.73e-3 * T .* (1 - exp(-1.263e-3 * T.^2)) .* ... % Specific heat capacity [J/kg K] 10.1051/0004-6361:20031746
    (1 + exp(-3 * T.^(1/2)) * 8.47e-3 .* T.^6 + ...
    2.0825e-7 * T.^4 .* exp(-4.97e-2 * T)));

% Thermal Expansivity [1/K]
a0        = 1.704e-4; % Ref thermal expansivity [1/K]
alpha_ice = a0 * T / BOD.Tice_0; % Thermal expansivity [1/K]

% Density
rho_ice   = (1-phi).*BOD.rhoIce_0.*(1-(T-BOD.Tice_0).*alpha_ice);  % Density [kg/m^3]


%%%%%%%%%%%%%%%%%%
% Calculate energy densities
%%%%%%%%%%%%%%%%%%
% Temperature energy density
U_T = rho_ice(1) * Cp_ice(1) * T(1) - rho_ice(2) * Cp_ice(2) * T(2);

% Latent heat energy density
U_m = BOD.rhoOcn*BOD.LH2O;

% Total
U = U_T + U_m;


%%%%%%%%%%%%%%%%%%
% Get radius
%%%%%%%%%%%%%%%%%%
rOcn_new = (dE/((4/3)*pi*U) + M.rOcn^3).^(1/3);



end










































