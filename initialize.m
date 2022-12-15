%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize Body
BOD = getBodyParameters(IN);

% Initialize Materials
[MAT,IN] = getMaterials(BOD,IN);

% Initialize Grid
[M, IN] = initializeGrid(IN,BOD);

% Initialize porosity
M.phi = zeros(size(M.r));

% Initialize Thermal Structure
[IN, M] = initializeThermal(IN,BOD,M,MAT);

% Initialize Time
[IN, M] = initializeTime(IN,M);

% Initialize Outputs
[OUT,MISC] = initializeOutputs(IN,M);