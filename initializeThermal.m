%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IN, M] = initializeThermal(IN,BOD,M,MAT)


%%%%%%%%%%%%%%%%%%%%%%%
% Create thermal structures
%%%%%%%%%%%%%%%%%%%%%%%
[M] = getInitialTemperatures(IN,BOD,M); % Initial surface temp;


%%%%%%%%%%%%%%%%%%%%%%%
% Initialize ocean parameters
%%%%%%%%%%%%%%%%%%%%%%%
M.vOcn  = (4/3)*pi*(M.rOcn^3-M.rSil^3); % Ocean volume
M.vH2O  = (4/3)*pi*(BOD.R^3-M.rSil^3);  % Hydrosphere volume
M.dm_dt = 0;  % Change in melt mass vs time

%%%%%%%%%%%%%%%%%%%%%%%
% Differentiation Properties
%%%%%%%%%%%%%%%%%%%%%%%
% Mass fraction rock differentiated on elements
M.fmDiff_H2O = zeros(1,M.Nz-1);
M.fmDiff_irn = ones(1,M.Nz-1);

% Leeched mass fraction (defined on elements!)
M.fmK = zeros(1,M.Nz-1);


%%%%%%%%%%%%%%%%%%%%%%%
% Material Arrays
%%%%%%%%%%%%%%%%%%%%%%%
% Track mass fractions of materials on elements
% [water, water melt, rock, rock melt, iron, iron melt]
M.mat.fm_s = zeros(6,M.Nz-1);
M.mat.iH2Osolid = 1; M.mat.iH2Omelts = 2;
M.mat.iSilSolid = 3; M.mat.iSilMelts = 4;
M.mat.iIrnSolid = 5; M.mat.iIrnMelts = 6;

M.mat.fm_s(M.mat.iH2Osolid,:) = IN.fm0_H2O;
M.mat.fm_s(M.mat.iSilSolid,:) = IN.fm0_sil;
M.mat.fm_s(M.mat.iIrnSolid,:) = IN.fm0_irn;

% Track volume fractions on elements
% [water, water melt, rock, rock melt, iron, iron melt]
% Need to initialize rho on nodes. Track for each
M.mat.rhoFull = zeros(6,M.Nz);
M.mat.rhoFull(M.mat.iH2Osolid,:) = MAT.H2O.s.rho0;
M.mat.rhoFull(M.mat.iH2Omelts,:) = MAT.H2O.m.rho0;
M.mat.rhoFull(M.mat.iSilSolid,:) = MAT.SIL.s.rho0;
M.mat.rhoFull(M.mat.iSilMelts,:) = MAT.SIL.m.rho0;
M.mat.rhoFull(M.mat.iIrnSolid,:) = MAT.IRN.s.rho0;
M.mat.rhoFull(M.mat.iIrnMelts,:) = MAT.IRN.m.rho0;
rhoFull_s = n2sVolumetric(M,M.mat.rhoFull);

% Bulk density
M.rho_s = sum(M.mat.fm_s./rhoFull_s).^-1;
M.rho   = s2nVolumetric(M,M.rho_s);

% Volume fractions
M.mat.fV_s = M.rho_s .* M.mat.fm_s ./ rhoFull_s; 


%%%%%%%%%%%%%%%%%%%%%%%
% Get Thermal Properties
%%%%%%%%%%%%%%%%%%%%%%%
M.k   = zeros(1,M.Nz);
M.rho = zeros(1,M.Nz);
M.Cp  = zeros(1,M.Nz);

[M] = getPressure(M,BOD,MAT);
[M] = getThermalProperties(M,IN,MAT);


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
























