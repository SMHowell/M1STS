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

if IN.radOn
    % Isotope Properties
    isotopeName = {'238U',     '235U',   '232Th',  '40K'  }; % Names of isotopes, purely for reference
    Piso        = [9.48e-5,    5.69e-4,  2.69e-5,  2.92e-5]; % Specific thermal power of isotope at present day [W/kg]
    tHalf       = [4.468,      0.7038,   14.05,    1.277  ]; % Half life [Gyr]
    fmIso       = [0.992745,   7.2e-3,   1.0,      1.17e-4]; % Abundance of isotope for each kg of parent element present [kg/kg]
    
    switch IN.primComp  % Elemental abundance (mass fraction) of parent element in rock [kg/kg]
        case 'CI'
            fm = [8e-9,   8e-9,   29e-9,  550e-6];
        case 'CM'
            fm = [12e-9,  12e-9,  41e-9,  370e-6];
        case 'CV'
            fm = [17e-9,  17e-9,  58e-9,  60e-6];
        case 'CO'
            fm = [18e-9,  18e-9,  80e-9,  360e-6];
        case 'CK'
            fm = [15e-9,  15e-9,  58e-9,  290e-6];
        case 'CR'
            fm = [13e-9,  13e-9,  42e-9,  315e-6];
    end
    
    % Calculate heating through time (more stable when left in Gyr)
    tPres = IN.age/IN.Gyr2s; % Present day [Ga]
    t     = M.t/IN.Gyr2s;    % Time of interest [Gyr]
    
    % Leeching fraction
    fmK_n = s2nMass(M,M.fmK); % Mass fraction leeched on nodes
    
    % Initial configuration for referencing
    fm0_sil = (1-(IN.fm0_H2O+IN.fm0_irn));
    
    % Specific power [W/kg]
    % For actual rocky interior
    PpRad_sil_1 = sum(Piso(1:3) .* fmIso(1:3) .* fm(1:3) .* exp(log(0.5)*(t-tPres)./tHalf(1:3)));
    PpRad_sil_2 = PpRad_sil_1 + (1-fmK_n) .* Piso(4) .* fmIso(4) .* fm(4) .* exp(log(0.5)*(t-tPres)./tHalf(4)); 
    PpRad_sil   = PpRad_sil_2 .* M.fm_sil./fm0_sil;    
    
    % Power [W] generated in each mantle layer
    m_sil     = M.fm_sil .* M.rho .* M.V;     % Mass of silicates in layer
    Prad_sil  = PpRad_sil .* m_sil / fm0_sil; % Radiogenic power 
    
    % Volumetric heat production [W/m^3]
    Hrad_sil  = Prad_sil./M.V; 
    
    % Ocean Heating ** NOTE! This is applied in getReservoirEnergy, not the
    % thermal solver **
    PpRad_ocn  = fmK_n .* Piso(4) .* fmIso(4) .* fm(4) .* exp(log(0.5)*(t-tPres)./tHalf(4)); % Heating rate [W/kg]
    PpRad_ocn  = PpRad_ocn .* M.fm_sil./fm0_sil;
    Prad_ocn   = PpRad_ocn .* m_sil / fm0_sil; % Radiogenic power 
    M.Prad_ocn = sum(Prad_ocn);  % Thermal power generation in the ocean from potassium
    
else
    M.Prad_ocn = 0;
    Hrad_sil   = zeros(size(M.H));
end

%%%%%%%%%%%%%%%%%%
% Tidal Heating
%%%%%%%%%%%%%%%%%%
if IN.tidalOn
    [M] = getIceRheology(M,BOD); % Get ice rheology
    
    % Get tidal heating rate on elements
    omega      = (2*pi)/BOD.tOrb; % Orbital forcing frequency
    Htidal_ice = (BOD.e0^2 * omega^2 * M.etaVE) ./ (1 + (omega^2 * M.etaVE.^2)/BOD.GH2O^2);
    
    % Set up array on elements
    H_s        = zeros(1,M.Nz-1);
    H_s(M.r_s>=M.rOcn & M.r_s<= BOD.R) = Htidal_ice(M.r_s>=M.rOcn & M.r_s<= BOD.R);
    
    % Heating rate on nodes [W/m^3]
    Htidal = s2nMass(M,H_s);
else
    Htidal = zeros(size(M.H));
end

%%%%%%%%%%%%%%%%%%
% Composite
%%%%%%%%%%%%%%%%%%
M.H = Hrad_sil+Htidal;
M.H(end) = 0;

end


























