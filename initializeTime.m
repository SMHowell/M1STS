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

% Initialize time array and timestep
M.t       = 0; % Time [s]
M.step    = 0; % Timestep
IN.outInd = 0; % Index for output to console
IN.outInd2 = 0; % Index for composition output to console

% Number of outputs
IN.NOut  = ceil(IN.tMax/IN.tOut)+1;
IN.NOut2 = ceil(IN.tMax/IN.tOut2)+1;

% Number of eruptions
IN.nErupt = 50;

end
























