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
% Differentiation
% Porosity
% Convection
% Silicate/iron melting
%%%%%

%%%%%%%%%%%%%%%%%%
% Model Timing
%%%%%%%%%%%%%%%%%%
IN.tMax  = '4.6 Gyr';   % Model run time from accretion, use 'yr' 'kyr' 'Myr' or 'Gyr'
IN.tOut  = '10 Myr';  % Output Frequency
IN.pltOn = 1;         % Plot output in real time
IN.movOn = 1;         % Save plot movie

%%%%%%%%%%%%%%%%%%
% Body Settings
%%%%%%%%%%%%%%%%%%
IN.body     = 'Europa';  % Target of interest
IN.lat      = [0, 90];   % Latitude (deg), can be a range e.g. [0, 90] to be averaged
IN.tidalOn  = 0;         % Turn on tidal heating
IN.radOn    = 1;         % Turn on radiogenic heating
IN.primComp = 'CI';      % Primordial Composition (CI, CM, CV, CO, CK, CR)
IN.age      = 4.6;       % Body age (present day) [Gyr]
IN.fm0_k    = 0.10;      % Max mass fraction of radioactive K leeched to water
IN.Tm_ocn   = 273;       % Melting temp of ocean [K] 

% Set only ONE of Hmax_x / fm0_x, and set the other to empty. If both are
% set, will default to Hmax_x setting. 
IN.Hmax_H2O = 130e3;     % Approx hydrosphere size after differentiation
IN.Hmax_irn = 500e3;     % Approx core size after differentiation

IN.fm0_H2O  = [];        % Mass fraction of water in primordial rock
IN.fm0_irn  = [];        % Mass fraction of iron in primordial rock

IN.T0_irn   = [];        % Initial isothermal core temperature [K]
                         % Set to empty to start at solar temperature. 
                         
%%%%%%%%%%%%%%%%%%
% Geometry
%%%%%%%%%%%%%%%%%%
IN.dz_H2O = 4e3;      % Minimum resolution in H2O layer [m]
IN.dz_sil = 50e3;     % Minimum resolution in rock [m]
IN.dz_irn = 50e3;     % Minimum resolution in core [m]

IN.H0_irn = 500e3;    % Initial iron core thickness   [m]

%%%%%%%%%%%%%%%%%%
% Run Model
%%%%%%%%%%%%%%%%%%
main;












