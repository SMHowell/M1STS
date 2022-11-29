%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = getHeating(BOD,M,IN)

%%%%%%%%%%%%%%%%%%
% Radiogenic
%%%%%%%%%%%%%%%%%%

% Calculate radiogenic heating, after Hussman et al. (2010) and references
% therein. BC and LYR contain model properties, M contains time,
% IN.primComp refers to the primordial composition, and fmK is a value between 0 and
% 1, where 1 represents total leeching to an ocean, and 0 represents no
% leeching at all.

% Isotope Properties 
isotopeName = {'238U',     '235U',   '232Th',  '40K'  }; % Names of isotopes, purely for reference
Piso        = [9.48e-5,    5.69e-4,  2.69e-5,  2.92e-5]; % Specific thermal power of isotope at present day [W/kg]
tHalf       = [4.468,      0.7038,   14.05,    1.277  ]; % Half life [Gyr]
fmIso       = [0.992745,   7.2e-3,   1.0,      1.17e-4]; % Abundance of isotope for each kg of parent element present [kg/kg]

switch IN.primComp  % Elemental abundance (mass fraction) of parent element in rock [kg/kg]
    case 'CI'
        fmSil = [8e-9,   8e-9,   29e-9,  (1-IN.fmK)*550e-6];
        fmOcn = [   0,      0,       0,      IN.fmK*550e-6];
    case 'CM'
        fmSil = [12e-9,  12e-9,  41e-9,  (1-IN.fmK)*370e-6];
        fmOcn = [   0,      0,       0,      IN.fmK*370e-6];
    case 'CV'
        fmSil = [17e-9,  17e-9,  58e-9,  (1-IN.fmK)*360e-6];
        fmOcn = [   0,      0,       0,      IN.fmK*360e-6];
    case 'CO'
        fmSil = [18e-9,  18e-9,  80e-9,  (1-IN.fmK)*360e-6];
        fmOcn = [   0,      0,       0,      IN.fmK*360e-6];
    case 'CK'
        fmSil = [15e-9,  15e-9,  58e-9,  (1-IN.fmK)*290e-6];
        fmOcn = [   0,      0,       0,      IN.fmK*290e-6];
    case 'CR'
        fmSil = [13e-9,  13e-9,  42e-9,  (1-IN.fmK)*315e-6];
        fmOcn = [   0,      0,       0,      IN.fmK*315e-6];
end

% Calculate heating through time (more stable when left in Gyr)
tPres = IN.age/IN.Gyr2s; % Present day [Ga]
t     = M.t/IN.Gyr2s;    % Time of interest [Gyr]

% Mantle Heating
PpRad_sil = sum(Piso .* fmIso .* fmSil .* exp(log(0.5)*(t-tPres)./tHalf)); % Heating rate [W/kg]
rho_s     = n2sVolumetric(M,M.rho);       % Density on elements
mRad      = sum(M.V_s(M.r_s<=M.rSil).*rho_s(M.r_s<=M.rSil)); % Mass of silicates and metals [kg]
Vsil      = (4/3)*pi*(M.rSil^3-M.rIrn^3); % Volume of silicates [m^3]
Hrad_sil  = PpRad_sil * mRad ./ Vsil;     % Volumetric heating rate in mantle [W/m^3]

% Ocean Heating ** NOTE! This is applied in getReservoirEnergy, not the
% thermal solver **
PpRad_ocn  = sum(Piso .* fmIso .* fmOcn .* exp(log(0.5)*(t-tPres)./tHalf)); % Heating rate [W/kg]
M.Prad_ocn = PpRad_ocn*mRad;  % Thermal power generation in the ocean from potassium

% Set up array on elements
H_s = zeros(1,M.Nz-1);
H_s(M.r_s>=M.rIrn & M.r_s<= M.rSil) = Hrad_sil;

% Heating rate on nodes [W/m^3]
Hrad_sil = s2nMass(M,H_s);



%%%%%%%%%%%%%%%%%%
% Tidal Heating
%%%%%%%%%%%%%%%%%%
[M] = getIceRheology(M,BOD); % Get ice rheology

% Get tidal heating rate on elements
omega      = (2*pi)/BOD.tOrb; % Orbital forcing frequency
Htidal_ice = (BOD.e0^2 * omega^2 * M.etaVE) ./ (1 + (omega^2 * M.etaVE.^2)/BOD.GH2O^2);

% Set up array on elements
H_s        = zeros(1,M.Nz-1);
H_s(M.r_s>=M.rOcn & M.r_s<= BOD.R) = Htidal_ice(M.r_s>=M.rOcn & M.r_s<= BOD.R);

% Heating rate on nodes [W/m^3]
H_tidal = s2nMass(M,H_s);



%%%%%%%%%%%%%%%%%%
% Composite
%%%%%%%%%%%%%%%%%%
M.H = Hrad_sil+H_tidal;


end


























