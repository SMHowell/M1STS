%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = thermalSolver(M,BOD,IN,MAT)

% Get tidal heating
M = getHeating(BOD,M,IN,MAT);

% Get surface heat flux
M = getSurfaceHeatFlux(M,BOD);

% Main thermal diffusion solve
refIndT = 2:M.Nz -1; % indices of solution
 
% Prefactor
A = - 2 * M.dt ./ (M.rho(refIndT).*M.Cp(refIndT).*M.r(refIndT).^2.*(M.dr(refIndT-1)+M.dr(refIndT)));

% Staggered node location
r_A = sqrt(M.r(refIndT).*M.r(refIndT-1));
r_B = sqrt(M.r(refIndT+1).*M.r(refIndT));

% Heat flux equality
k_A = (M.k(refIndT).*M.k(refIndT-1).*(M.r(refIndT)-M.r(refIndT-1)))./(M.k(refIndT).*(r_A-M.r(refIndT-1))+M.k(refIndT-1).*(M.r(refIndT)-r_A));
k_B = (M.k(refIndT+1).*M.k(refIndT).*(M.r(refIndT+1)-M.r(refIndT)))./(M.k(refIndT+1).*(r_B-M.r(refIndT))+M.k(refIndT).*(M.r(refIndT+1)-r_B));

% Set thermal conductivity on ocean and reservoir boundaries
[k_A, k_B] = boundaryConductivity(M,k_A,k_B);

% Staggered node heat fluxes 
q_A = - M.Nu(refIndT-1) .* k_A .* (M.T(refIndT)-M.T(refIndT-1)) ./ (M.dr(refIndT-1));
q_B = - M.Nu(refIndT)   .* k_B .* (M.T(refIndT+1)-M.T(refIndT)) ./ (M.dr(refIndT));

% Staggered node heat flows
Q_A = r_A.^2 .* q_A;
Q_B = r_B.^2 .* q_B;

% Apply heat flow boundary and conditions
Q_B(end) = M.qLoss * r_B(end)^2;

% Temperature change
M.dT(refIndT) = A .* (Q_B - Q_A) + (M.dt./(M.rho(refIndT) .* M.Cp(refIndT))) .* M.H(refIndT);

% Update temperatures
M.T      = M.T + M.dT;

% Apply boundary conditions to ghost nodes
M.Tsurf  = M.T(end-1); % Track new surface temperature
M.T(end) = M.Tsurf;    % Zero out dT in space
M.T(1)   = M.T(2);     % dT/dr = 0 at r = 0

if any(isnan(M.T))
    imBroken
end

end























