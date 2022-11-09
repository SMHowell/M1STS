%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% reservoir freezing and evolution.
% Sam Howell, Elodie Lesage
% samuel.m.howell@jpl.nasa.gov
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT,MISC] = initializeOutputs(IN,M)

% Initialize and add the arrays to track here
OUT.t     = zeros(1,IN.NOut); % Time of output
OUT.dm_dt = zeros(1,IN.NOut); % Freezing rate [kg/s]
OUT.Tsurf = zeros(1,IN.NOut); % Surface Temperature
OUT.Dice  = zeros(1,IN.NOut); % Ice thickness

OUT.vRes  = zeros(1,IN.NOut); % Reservoir volume
OUT.rRes  = zeros(1,IN.NOut); % Reservoir radius
OUT.zRes  = zeros(1,IN.NOut); % Reservoir Depth
OUT.fV    = zeros(1,IN.NOut); % Reservoir frozen fraction
OUT.Tmelt = zeros(1,IN.NOut); % Reservoir melting temperature

MISC      = []; % Miscellaneous array collector

end