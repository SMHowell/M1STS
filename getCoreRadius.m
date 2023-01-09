%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,rIoc_new] = getCoreRadius(M,BOD,IN,dE,MAT)

%%%%%%%%%%%%%%%%%%
% Iron Thermal Properties
%%%%%%%%%%%%%%%%%%
% Get props at appropraite temps
% Handle degenerate case where interface is already at melting temp
if M.T(M.iIocTop+1) < MAT.IRN.Tm0
    T_avg = (MAT.IRN.Tm0-M.T(M.iIocTop+1))/log(MAT.IRN.Tm0/M.T(M.iIocTop+1));
else
    T_avg = MAT.IRN.Tm0;
end
T   = [MAT.IRN.Tm0, T_avg];

% Heat capacity
Cp_irn = [MAT.IRN.s.Cp0,MAT.IRN.s.Cp0];

% Density
rho_irn   = [MAT.IRN.s.rho0,MAT.IRN.s.rho0];  % Density [kg/m^3]


%%%%%%%%%%%%%%%%%%
% Calculate energy densities
%%%%%%%%%%%%%%%%%%
% Temperature energy density
U_T = rho_irn(1) * Cp_irn(1) * T(1) - rho_irn(2) * Cp_irn(2) * T(2);

% Latent heat energy density
U_m = MAT.IRN.m.rho0*MAT.IRN.L;

% Total
U = U_T + U_m;


%%%%%%%%%%%%%%%%%%
% Get radius
%%%%%%%%%%%%%%%%%%
rIoc_new = (-dE/((4/3)*pi*U) + M.rIoc^3).^(1/3);



end










































