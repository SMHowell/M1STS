%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% reservoir freezing and evolution.
% Sam Howell, Elodie Lesage
% samuel.m.howell@jpl.nasa.gov
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = initializeGrid(IN,BOD)

% Define useful values
rBot     = BOD.R - 2*IN.H0; % Radius at the base of the model (2 ice shell depths) [m]
rIce     = BOD.R - IN.H0;   % Radius at the base of the ice shell [m]
M.refInd = 1:IN.Nz-1;       % A useful index for referencing values 
M.rOcnTop = rIce;           % Depth to ocean surface


% Made radial grid extending to ghost node in space
dr0  = (2*IN.H0/(IN.Nz-1)); % Approximate dr
M.r  = [linspace(rBot,rIce,floor((IN.Nz-1)/2)),linspace(rIce+dr0,BOD.R,ceil((IN.Nz-1)/2)), BOD.R+dr0]; % Radius from center of body [m]
M.dr = diff(M.r); % Element size [m]
M.z  = BOD.R - M.r; % Depth [z]

% Define the location of elements in thermal calculations. Note that M.r_s
% converges to the midpoint between nodes when ice thickness << body radius
M.r_s = sqrt(M.r(M.refInd) .* M.r(M.refInd+1)); 
M.z_s = BOD.R - M.r_s; % Depth [z]

% Define element (staggered node) volumes, surface areas, and length scales
M.V_s = (4/3) * pi * (M.r(M.refInd+1).^3-M.r(M.refInd).^3);
M.A   = 4 * pi * M.r.^2;
M.L   = M.V_s./M.A(M.refInd+1);  % characteristic length scale

end