%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OUT,MISC] = initializeOutputs(IN,M)

% Initialize and add the arrays to track here
OUT.t     = zeros(1,IN.NOut); % Time of output
OUT.dm_dt = zeros(1,IN.NOut); % Freezing rate [kg/s]
OUT.Tsurf = zeros(1,IN.NOut); % Surface Temperature
OUT.Dice  = zeros(1,IN.NOut); % Ice thickness
OUT.DH2O  = zeros(1,IN.NOut); % H2O thickness

% Track start and end of differentiation
OUT.diffH2Ostart = [];
OUT.diffH2Ostop  = [];
OUT.diffIrnStart = [];
OUT.diffIrnStop  = [];


MISC      = []; % Miscellaneous array collector

end