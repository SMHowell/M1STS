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
        BOD.A       = 0.55;    % bolometric albedo
        BOD.eps     = 0.9;     % Emissivity
        BOD.ob      = 3;       % obliquity [deg]
        BOD.RParent = 7.78e11; % orbital radius of system from sun [m]
        BOD.tOrb    = 306800;  % orbital period
        BOD.e0      = 1.5e-5;  % Tidal strain amplitude
        
    otherwise
        error('Error. %s is not definited in getBodyParameters.m.',IN.body)
        
end

% Derivitive parameters
BOD.V       = (4/3)*pi*BOD.R^3; % Body volume [m^3]
BOD.rhoBulk = BOD.m/BOD.V;      % Bulk density [kg/m^3] 

%%%%%%%%%%%%%%%%%%%%%%%
% Parameters and Materials common to all bodies
%%%%%%%%%%%%%%%%%%%%%%%
% General Properties
BOD.vfmCr      = 0.40;   % Melt fraction required for disaggregation (~0.4 - 0.6):
                         % https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2002JE001943
BOD.aMelt      = 30;     % Exponent on melt-fraction dependence of viscosity: 
                         % soest.hawaii.edu/GG/FACULTY/smithkonter/GG631/other/HirthKohlstedt_2000.pdf
                        
% H2O Properties
BOD.GH2O       = 3e9;    % Elastic modulus of ice (GPa)
BOD.EaH2O      = 59.4e3; % Activation Energy [J mol-1] 
BOD.LH2O       = 330e3;  % Latent heat of fusion [kJ/kg]
BOD.etaOcn_0   = 1e-3;   % Viscosity of melt water [Pa s]

BOD.rhoOcn     = 1000;   % Density [kg/m^3]
BOD.CpOcn      = 4184;   % Specific heat capacity [J/kg K]
BOD.rhoIce_0   = 917;    % Ice reference density [kg/m^3]
BOD.CpIce_0    = 2107;   % Specific heat capacity [J/kg K]
BOD.Tice_0     = 273;    % Reference temeprature for laws [K]
 
% Rock: https://www.sciencedirect.com/science/article/pii/S0031920113000289
% using forsterite (Mg-olivine) as representative
BOD.TmSil      = 1300;  % Melting temp of rock [K] 
BOD.kSil       = 4.0;   % Thermal conductivity [W/mK]
BOD.Lsil       = 418e3; % Latent heat of fusion [kJ/kg]
BOD.rhoSil_0   = 3275;  % Reference Density [kg/m^3]
BOD.Tsil_0     = 300;   % Reference temeprature for laws [K] 
                         
BOD.TmIrn      = 1260;  % Melting temp of iron [K] 
BOD.kIrn_0     = 27.0;  % Thermal conductivity [W/mK]
BOD.rhoIrn_0   = 5428;  % Density [kg/m^3]
BOD.rhoIrn_m_0 = 5240;  % Density [kg/m^3]
BOD.CpIrn_0    = 450;   % Specific heat capacity [J/kg K]
BOD.CpIrn_m_0  = 820;   % Specific heat capacity of melt [J/kg K]
BOD.Lirn       = 247e3; % Latent heat of fusion [kJ/kg]


end