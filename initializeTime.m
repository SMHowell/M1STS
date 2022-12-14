%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IN,M] = initializeTime(IN,M)

% Book keeping on inputs
IN = parseTime(IN);

% Book keep timing
M.simuTime_old = 0; % Simulation time
M.realTime_old = 0; % Elapsed time
tic; 

% Initialize time array and timestep
M.t       = 0; % Time [s]
M.step    = 0; % Timestep
IN.outInd = 0; % Index for output to console

% Number of outputs
IN.NOut = ceil(IN.tMax/IN.tOut)+1;

end
























