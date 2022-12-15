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

%%%%%%%%%%%%%%%%%%%%%%%
% Body Properties
%%%%%%%%%%%%%%%%%%%%%%%
% Add bodies as new cases below
switch IN.body
    case 'Europa'
        BOD.g       = 1.315;   % Gravity [m/s^2]
        BOD.R       = 1561e3;  % Body radius [m]
        BOD.m       = 4.80e22; % Body mass [kg]
        BOD.BoloA   = 0.55;    % bolometric albedo
        BOD.eps     = 0.9;     % Emissivity
        BOD.ob      = 3;       % obliquity [deg]
        BOD.RParent = 7.78e11; % orbital radius of system from sun [m]
        BOD.tOrb    = 306800;  % orbital period
        BOD.e0      = 1.5e-5;  % Tidal strain amplitude
        
    otherwise
        error('Error. %s is not definited in getBodyParameters.m.',IN.body)
        
end

% Derivitive parameters
BOD.Asurf   = 4*pi*BOD.R^2;     % Body surface area [m^3]
BOD.V       = (4/3)*pi*BOD.R^3; % Body volume [m^3]
BOD.rhoBulk = BOD.m/BOD.V;      % Bulk density [kg/m^3] 
BOD.omega   = (2*pi)/BOD.tOrb;  % Orbital forcing frequency [rad/s]

end