%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = getSurfaceHeatFlux(M,BOD)

% Current surface temp 
Tsurf = M.Tsurf; 

% Cosntants
sigma  = 5.67e-8;     % Stefan Boltzamnn constant [W m^-2 K^-4]
m      = 18.01528e-3; % Molar weight of H2O [kg mol^-1]
R      = 8.314;       % Ideal gas constant [J mol^-1 K^-1]

%%%%%%%%%%%%%%%%%%%%%%%
% Emissive heat flux
%%%%%%%%%%%%%%%%%%%%%%%
M.qEmis = BOD.eps * sigma * Tsurf.^4;

%%%%%%%%%%%%%%%%%%%%%%%
% Sublimation heat flux
%%%%%%%%%%%%%%%%%%%%%%%
% Calculate latent heat of sublimation
% Polynomial curve fits to Table 2.1. R. R. Rogers; M. K. Yau (1989).
% BOD.BoloA Short Course in Cloud Physics (3rd ed.). Pergamon Press. p. 16.
% ISBN 0-7506-3215-1.
L = (2834.1 - 0.29*Tsurf - 0.004*Tsurf.^2)*1e3;

% Calculate the partial pressure of ice
% https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/93GL00105
% NOTE log = log10 in this paper, not ln
pA = -2663.5; % empirical coefficient
pB = 12.537;  % empirical coefficient
P = 10.^(pA./Tsurf + pB); % partial pressure [Pa]

% heat flux
M.qSub = L .* P .* sqrt(m./(2 * pi * R * Tsurf));

%%%%%%%%%%%%%%%%%%%%%%%
% Required Geologic Heat Flux [W/m^2]
%%%%%%%%%%%%%%%%%%%%%%%
qIn   = M.qSol;           % Paid into heat budget
qOut  = M.qEmis + M.qSub; % Withdrawn from heat budget
M.qLoss = qOut-qIn;       % Anomalous heat flux

end








