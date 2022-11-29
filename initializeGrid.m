%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = initializeGrid(IN,BOD)
%%%%%%%%%%%%%%%%%%%%%%
% Model Grid
%%%%%%%%%%%%%%%%%%%%%%
M.rOcn = BOD.R - IN.H0_ice;  % Outer radius of ocean surface
M.rSil = BOD.R - IN.H0_H2O;  % Outer radius of seafloor
M.rIrn = IN.H0_irn;          % Outer radius of core

r_ice = linspace(M.rOcn, BOD.R,  ceil((BOD.R -M.rOcn)/IN.dz_H2O));
r_ocn = linspace(M.rSil, M.rOcn, ceil((M.rOcn-M.rSil)/IN.dz_H2O));
r_sil = linspace(M.rIrn, M.rSil, ceil((M.rSil-M.rIrn)/IN.dz_sil));
r_irn = linspace(0,      M.rIrn, ceil(M.rIrn/IN.dz_irn));

% Radial array with space ghost node 1 km off the surface
M.r   = sort(unique([r_ice,r_ocn,r_sil,r_irn,BOD.R+IN.dz_H2O]));
M.dr  = diff(M.r);             % Element height [m]
M.z   = BOD.R - M.r;           % Depth [m]

M.Nz = numel(M.r);            % Number of nodes in vertical direction
M.ind = 1:M.Nz-1;             % A useful index for referencing values 

% Define the location of elements in thermal calculations. Note that M.r_s
% converges to the midpoint between nodes when ice thickness << body radius
M.r_s = sqrt(M.r(M.ind) .* M.r(M.ind+1));  % Radius from center of body [m]
M.z_s = BOD.R - M.r_s;                     % Depth [m]

% Define element (staggered node) volumes, surface areas, and length scales
M.V_s = (4/3) * pi * (M.r(M.ind+1).^3-M.r(M.ind).^3);
M.A   = 4 * pi * M.r.^2;
M.L   = M.V_s./M.A(M.ind+1);  % characteristic length scale

%%%%%%%%%%%%%%%%%%%%%%
% The volumetric weightings are pre-computed to save time in later
% interpolations
%%%%%%%%%%%%%%%%%%%%%%
M.Vf1 = (M.r_s(M.ind).^3 - M.r(M.ind).^3)./  (M.r(M.ind+1).^3 - M.r(M.ind).^3); % Weighting for element i
M.Vf2 = (M.r(M.ind+1).^3 - M.r_s(M.ind).^3)./(M.r(M.ind+1).^3 - M.r(M.ind).^3); % Weighting for element i+1

end














