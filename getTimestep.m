%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dt] = getTimestep(M)
% Courant-Friedreichs-Lewy condition
CFL = 1/3;

% Thermal diffusion timestep
% Estimate K on staggered nodes (imperfectly)
[K_s] = n2sVolumetric(M,M.K);
dt_th = min(M.L.^2 ./ K_s);

% Radiative emission timestep (approx)
dt_emis = (M.rho(end-1) * M.Cp(end-1) * M.Tsurf * M.L(end-1))/M.qSol;

% Update Timestep
dt = CFL * min(dt_th,dt_emis);

end