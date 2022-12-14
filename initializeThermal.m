%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IN, M] = initializeThermal(IN,BOD,M)


%%%%%%%%%%%%%%%%%%%%%%%
% Create thermal structures
%%%%%%%%%%%%%%%%%%%%%%%
[M] = getInitialTemperatures(IN,BOD,M); % Initial surface temp;


%%%%%%%%%%%%%%%%%%%%%%%
% Initialize melts
%%%%%%%%%%%%%%%%%%%%%%%
% Melt fractions for continuous layers (e.g. ocean, outer core)
% *volume* fraction (on elements)
M.vfm = zeros(1,M.Nz-1);         
M.vfm(M.r>=M.rSil & M.r< M.rOcn) = 1; % Impose ocean
M.dm_dt = 0;  % Change in melt mass vs time

% Partial melt fractions
% *volume* fraction (on elements)
M.vpfm  = zeros(1,M.Nz-1);         


%%%%%%%%%%%%%%%%%%%%%%%
% Initialize reservoir/ocean parameters
%%%%%%%%%%%%%%%%%%%%%%%
M.vOcn    = (4/3)*pi*(M.rOcn^3-M.rSil^3); % Ocean volume
M.vH2O    = (4/3)*pi*(BOD.R^3-M.rSil^3);  % Hydrosphere volume


%%%%%%%%%%%%%%%%%%%%%%%
% Differentiation Properties
%%%%%%%%%%%%%%%%%%%%%%%
% Mass fraction rock differentiated on elements
M.fmDiff_H2O = zeros(1,M.Nz-1);
M.fmDiff_irn = ones(1,M.Nz-1);

% Leeched mass fraction (defined on elements!)
M.fmK = zeros(1,M.Nz-1);


%%%%%%%%%%%%%%%%%%%%%%%
% Get Thermal Properties
%%%%%%%%%%%%%%%%%%%%%%%
M.k   = zeros(1,M.Nz);
M.rho = zeros(1,M.Nz);
M.Cp  = zeros(1,M.Nz);

[M] = getPressure(M,BOD);
[M] = getThermalProperties(M,BOD,IN);


%%%%%%%%%%%%%%%%%%%%%%%
% Nusselt Number on Elements
%%%%%%%%%%%%%%%%%%%%%%%
M.Ra_cr = 0; % Critical Rayleigh Number
M.Nu    = ones(1,M.Nz-1);

%%%%%%%%%%%%%%%%%%%%%%%
% Tidal heating on nodes
%%%%%%%%%%%%%%%%%%%%%%%
M.H = zeros(1,M.Nz);


end
























