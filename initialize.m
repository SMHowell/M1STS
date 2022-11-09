%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% reservoir freezing and evolution.
% Sam Howell, Elodie Lesage
% samuel.m.howell@jpl.nasa.gov
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Body
BOD = getBodyParameters(IN);

% Initialize Grid
M = initializeGrid(IN,BOD);

% Initialize porosity
M.phi = zeros(size(M.r));

% Initialize Thermal Structure
M = initializeTemp(IN,BOD,M);

% Initialize Time
[IN, M] = initializeTime(IN,M);

% Initialize Outputs
[OUT,MISC] = initializeOutputs(IN,M);