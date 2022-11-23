%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = getPressure(M,BOD)

 % Guestimate -- *** REPLACE THIS ***
rhoOcn   = 1000;  % Density [kg/m^3]
rho0_sil = 3275;  % Reference Density [kg/m^3]

P_min = (BOD.R-M.rSil)  * rhoOcn * BOD.g;
P_max = (M.rSil-M.rIrn) * rho0_sil * BOD.g;
M.P   = linspace(P_min,P_max,M.Nz);

end









































