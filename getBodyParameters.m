%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% reservoir freezing and evolution.
% Sam Howell, Elodie Lesage
% samuel.m.howell@jpl.nasa.gov
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
                
    otherwise
        error('Error. %s is not definited in getBodyParameters.m.',IN.body)
        
end


end