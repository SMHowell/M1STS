%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IN = parseTime(IN)

% seconds in time unit
IN.yr2s  = 365.25*24*3600;  % Handy dandy conversion
IN.kyr2s = 1e3 * IN.yr2s;   % Handy dandy conversion
IN.Myr2s = 1e6 * IN.yr2s;   % Handy dandy conversion
IN.Gyr2s = 1e9 * IN.yr2s;   % Handy dandy conversion

% Parse final time
timeStr = strsplit(IN.tMax);
val     = str2double(timeStr{1});
unit    = timeStr{2};

switch unit
    case 'yr'
        conversion = IN.yr2s;
    case 'kyr'
        conversion = IN.kyr2s;
    case 'Myr'
        conversion = IN.Myr2s;
    case 'Gyr'
        conversion = IN.Gyr2s;
    otherwise
        error('Incorrect timescale selected.');
end
IN.tMax = val*conversion;

% Parse output time
timeStr = strsplit(IN.tOut);
val     = str2double(timeStr{1});
unit    = timeStr{2};

switch unit
    case 'yr'
        conversion = IN.yr2s;
    case 'kyr'
        conversion = IN.kyr2s;
    case 'Myr'
        conversion = IN.Myr2s;
    case 'Gyr'
        conversion = IN.Gyr2s;
    otherwise
        error('Incorrect timescale selected.');
end
IN.tOut = val*conversion;

% Get number of output timesteps
IN.outN = ceil(IN.tMax / IN.tOut)+1; 


% Parse emplacement delay
timeStr = strsplit(IN.tRes);
val     = str2double(timeStr{1});
unit    = timeStr{2};

switch unit
    case 'yr'
        conversion = IN.yr2s;
    case 'kyr'
        conversion = IN.kyr2s;
    case 'Myr'
        conversion = IN.Myr2s;
    case 'Gyr'
        conversion = IN.Gyr2s;
    otherwise
        error('Incorrect timescale selected.');
end
IN.tRes = val*conversion;



end























