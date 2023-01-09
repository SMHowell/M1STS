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
IN.tMax  = '1 Myr';   % Model run time, use 'yr' 'kyr' 'Myr' or 'Gyr'
IN.tOut  = '10 kyr';  % Output Frequency
IN.tOut2 = '0.1 kyr';  % Output Frequency 2
IN.pltOn = 1;         % Plot output in real time
IN.movOn = 0;         % Save plot movie

%%%%%%%%%%%%%%%%%%
% Body Settings
%%%%%%%%%%%%%%%%%%
IN.body = 'Europa';   % Target of interest
IN.lat  = 0;          % Latitude of simulation (deg)

%%%%%%%%%%%%%%%%%%
% Geometry
%%%%%%%%%%%%%%%%%%
IN.H0  = 10e3;        % Initial ice shell thickness [m]

IN.Nz  = 200;         % Number of nodes in vertical direction
IN.Nx  = 100;         % Number of nodes in hotizontal direction (not implemented)

% NOTE: You always want at least 3 full elements between the surface and
% top of your reservoir, and b/w the bottom of the reservoir and ice-ocean
% interface because of thermal conductivity smoothing for convective water
% bodies. It may break when the number of reservoir elements is <= 3.
% Maybe not.
IN.rRes    = 2e3;         % Initial reservoir radius [m]
IN.zResTop = 2e3;         % Initial reservoir top depth [m]
IN.tRes    = '10 kyr';   % Time after which to start emplacement


%%%%%%%%%%%%%%%%%%
% Thermal Properties
%%%%%%%%%%%%%%%%%%
% IN.Tm_ocn = 273;   % Melting temp of ocean [K] - NOTE: This should set 
% equal to the initial reservoir melting temperature 
IN.CpOcn  = 4184;  % Specific heat capacity [J/kg K]
IN.L      = 330e3; % Latent heat of fusion [kJ/kg]
IN.X      = 5e-10; % Liquid water compressibility [Pa^-1]

%%%%%%%%%%%%%%%%%%
% Composition
%%%%%%%%%%%%%%%%%%
IN.simu   = 14;      % simulation number (see CompData file)
IN.Ncomp  = 8;       % Number of species tracked

% Number of eruptions
IN.nErupt = 1e3;

%%%%%%%%%%%%%%%%%%
% Run Model
%%%%%%%%%%%%%%%%%%
main;












