%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

%%%%% TBD:
% Horizontal solver
% Porosity
% Tidal heat gen
% Convection
%%%%%

%%%%%%%%%%%%%%%%%%
% Model Timing
%%%%%%%%%%%%%%%%%%
IN.tMax  = '4 Gyr';   % Model run time, use 'yr' 'kyr' 'Myr' or 'Gyr'
IN.tOut  = '0.1 Myr';   % Output Frequency
IN.pltOn = 1;         % Plot output in real time
IN.movOn = 0;         % Save plot movie


%%%%%%%%%%%%%%%%%%
% Body Settings
%%%%%%%%%%%%%%%%%%
IN.body     = 'Europa';  % Target of interest
IN.lat      = [0, 90];   % Latitude (deg), can be a range e.g. [0, 90]
IN.primComp = 'CI';      % Primordial Composition (CI, CM, CV, CO, CK, CR)
IN.age      = 4.5;       % Body age [Gyr]
IN.fmK      = 0.05;      % Mass fraction of radioactive K leeched to ocean

%%%%%%%%%%%%%%%%%%
% Geometry
%%%%%%%%%%%%%%%%%%
IN.dz_H2O = 1e3;      % Minimum resolution in H2O layer [m]
IN.dz_sil = 50e3;     % Minimum resolution in rock [m]
IN.dz_irn = 20e3;     % Minimum resolution in core [m]

IN.H0_ice = 130e3;    % Initial ice shell thickness   [m]
IN.H0_H2O = 130e3;    % Initial hydrosphere thickness [m]
IN.H0_irn = 500e3;    % Initial iron core thickness   [m]

IN.T0_irn = 1260;     % Initial isothermal core temperature [K]
IN.Tm_ocn = 273;      % Melting temp of ocean [K] 

% NOTE: You always want at least 3 full elements between the surface and
% top of your reservoir, and b/w the bottom of the reservoir and ice-ocean
% interface because of thermal conductivity smoothing for convective water
% bodies. It may break when the number of reservoir elements is <= 3.
% Maybe not.
IN.resOn   = 0;          % Enable Ice shell reservoir
IN.rRes    = 5e3;        % Initial reservoir radius [m]
IN.zResTop = 10e3;       % Initial reservoir top depth [m]
IN.tRes    = '100 kyr';  % Time after which to start emplacement


%%%%%%%%%%%%%%%%%%
% Run Model
%%%%%%%%%%%%%%%%%%
main;












