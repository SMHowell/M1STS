%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% reservoir freezing and evolution.
% Sam Howell, Elodie Lesage
% samuel.m.howell@jpl.nasa.gov
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IN,M] = initializeTime(IN,M)

% Book keeping on inputs
IN = parseTime(IN);

% Initialize time array and timestep
M.t       = 0; % Time [s]
M.step    = 0; % Timestep
IN.outInd = 0; % Index for output to console

% Number of outputs
IN.NOut = ceil(IN.tMax/IN.tOut)+1;

end
























