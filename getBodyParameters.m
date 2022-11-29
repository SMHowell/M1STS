%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BOD] = getBodyParameters(IN)

% Add bodies as new cases below
switch IN.body
    case 'Europa'
        BOD.g       = 1.315;   % Gravity [m/s^2]
        BOD.R       = 1561e3;  % Body radius [m]
        BOD.m       = 4.80e22; % Body mass [kg]
        BOD.A       = 0.55;    % bolometric albedo
        BOD.eps     = 0.9;     % Emissivity
        BOD.ob      = 3;       % obliquity [deg]
        BOD.RParent = 7.78e11; % orbital radius of system from sun [m]
        BOD.tOrb    = 306800;  % orbital period
        BOD.e0      = 1.5e-5;  % Tidal strain amplitude
        
    otherwise
        error('Error. %s is not definited in getBodyParameters.m.',IN.body)
        
end

% Parameters common to all bodies
% Ice Properties
BOD.GH2O  = 3e9;    % Elastic modulus of ice (GPa)
BOD.EaH2O = 59.4e3; % Activation Energy [J mol-1] 

end