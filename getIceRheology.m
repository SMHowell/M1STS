%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric thermal solver
% Conservative finite differences, explicit forward in time
% Samuel.M.Howell@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = getIceRheology(M,BOD)

% Constants
R = 8.314;  % Ideal gas constant [J mol-1 K-1]
omega = (2*pi)/BOD.tOrb; % Orbital forcing frequency

% Material Inputs 
T0   = 273;  % Reference temperature for eta0 [K]
eta0 = BOD.GH2O/omega; % Reference viscosity at T0 [Pa s]

% Get creep viscosity
etaDiff = eta0 * exp((BOD.EaH2O/(R*T0))*(T0./M.T-1));

% Get visco-elastic effective viscosity
Z       = BOD.GH2O*M.dt./(BOD.GH2O*M.dt+etaDiff); % Elasticity correction
M.etaVE = Z.*etaDiff;

end