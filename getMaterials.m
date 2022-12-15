%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MAT, IN] = getMaterials(BOD,IN)

%%%%%%%%%%%%%%%%%%%%%%%
% General Properties
%%%%%%%%%%%%%%%%%%%%%%%
MAT.R = 8.314;  % Ideal gas constant [J mol-1 K-1]

%%%%%
% H2O
%%%%%
% General
MAT.H2O.Tm0 = IN.Tm_ocn; % Melting temperature [K]
MAT.H2O.L   = 330e3;  % Latent heat of fusion [kJ/kg]

% Solids
MAT.H2O.s.nu    = 0.33;   % Poisson's ratio
MAT.H2O.s.Emod  = 9e9;    % Young's modulus [Pa]
MAT.H2O.s.Gmod  = MAT.H2O.s.Emod/(2*(1+MAT.H2O.s.nu));   % Elastic shear modulus [Pa]
MAT.H2O.s.Kmod  = MAT.H2O.s.Emod/(3*(1-2*MAT.H2O.s.nu)); % Bulk modulus [Pa]
MAT.H2O.s.Ea    = 59.4e3; % Activation Energy [J mol-1]
MAT.H2O.s.T0    = 273;    % Reference temperature for props [K]
MAT.H2O.s.eta0  = MAT.H2O.s.Gmod/BOD.omega; % Reference viscosity at T0 [Pa s]
MAT.H2O.s.rho0  = 917;    % Reference density [kg/m^3]
MAT.H2O.s.Cp0   = 2107;   % Specific heat capacity [J/kg K]

% Melts
MAT.H2O.m.Kmod  = 2.2e9;  % Bulk modulus [Pa]
MAT.H2O.m.T0    = 273;    % Reference temperature for props [K]
MAT.H2O.m.eta0  = 1e-3;   % Viscosity of melt  [Pa s]
MAT.H2O.m.rho0  = 1000;   % Density [kg/m^3]
MAT.H2O.m.Cp0   = 4184;   % Specific heat capacity [J/kg K]




%%%%%
% Rock
%%%%%
switch IN.silComp
    case 'Fo'
        %%%%%
        % Forsterite (Mg-Olivine)
        %%%%%
        % Mech. https://www.mdpi.com/2075-163X/9/12/787
        % Therm. https://www.sciencedirect.com/science/article/pii/S0031920113000289
        % Rheo. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013JB010473
        % Melting. https://www.sciencedirect.com/science/article/pii/0031920181900844
        % Magma. https://www.geochemsoc.org/files/4414/1258/4210/SP-1_047-058_Herzberg.pdf
        % http://www2.tulane.edu/~sanelson/Natural_Disasters/volcan&magma.htm#:~:text=Thus%2C%20basaltic%20magmas%20tend%20to,times%20more%20viscous%20than%20water.
        % Scratch
        Emod = [274.4, 153.2, 170.9]*1e9; % Triaxial Young's modulus [Pa]
        nu   = [0.23, 0.20, 0.13, 0.24, 0.14, 0.29]; % Poisson's ratio tensor
        A0   = 1.042e17; % [Pa s / m^2] Diffusion prefactor
        d    = 1e-6;   % Grain size [m]
        p    = 2;      % Grain size exponent
        P0     = 100e6;   % Reference pressure [Pa]
        
        % General
        MAT.SIL.L     = 506e3;  % Latent heat of fusion [kJ/kg]
        MAT.SIL.Tm0   = 2163;   % Reference melting temperature at surface [K]
        % 2163 K for forsterite
        % 1270 - 1470 K for basalt (low viscosity)
        % 930 - 1070 K for Rhyolite (high viscosity)
        
        % Solids
        MAT.SIL.s.nu    = mean(nu);     % Poisson's ratio
        MAT.SIL.s.Emod  = mean(Emod);   % Young's modulus [Pa]
        MAT.SIL.s.Gmod  = MAT.SIL.s.Emod/(2*(1+MAT.SIL.s.nu));   % Elastic shear modulus [Pa]
        MAT.SIL.s.Kmod  = MAT.SIL.s.Emod/(3*(1-2*MAT.SIL.s.nu)); % Bulk modulus [Pa]
        MAT.SIL.s.T0    = 2163;   % Reference temperature for props [K]
        MAT.SIL.s.Ea    = 469e3;  % Activation Energy [J mol-1]
        MAT.SIL.s.Va    = 7.5e-6; % Activation Volume [m^3 mol-1]
        MAT.SIL.s.eta0  = A0*d^p*exp((MAT.SIL.s.Ea + P0 * MAT.SIL.s.Va)./(MAT.R*MAT.SIL.s.T0)); % Reference viscosity at T0 [Pa s]
        MAT.SIL.s.k0    = 2.47;   % Referemce thermal conductivity [W m-1 K-1]
        MAT.SIL.s.rho0  = 3275;   % Reference density [kg/m^3]
        MAT.SIL.s.Cp0   = 875.5;  % Specific heat capacity [J/kg K]
        
        % Melts
        MAT.SIL.m.Kmod  = 22.6e9; % Bulk modulus [Pa]
        MAT.SIL.m.eta0  = 1e3;    % Viscosity of melt [Pa s]
        MAT.SIL.m.rho0  = 2520;   % Density [kg/m^3]
        MAT.SIL.m.Cp0   = 1600;   % Specific heat capacity [J/kg K]
        
    case 'Fa'
        %%%%%
        % Fayalite (Fe-Olivine)
        %%%%%
        % https://link.springer.com/article/10.1007/BF00203203
        % https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JZ072i016p04235
        % Rheo. file:///Users/smhowell/Desktop/1-s2.0-0012821X96001549-main.pdf
        % Density. http://www.webmineral.com/data/Fayalite.shtml
        % Expansion. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC117344/
        % https://www.osti.gov/biblio/6720778
        % Scratch
        Kmod = 127.9e9;
        Gmod = 50.3e9;
        Emod = 9*Gmod*Kmod/(Gmod+3*Kmod);
        nu   = Emod/(2*Gmod) -1;
        
        n      = 3.5;     % Stress exponent
        d      = 1e-6;    % Grain size [m]
        A0_dry = 2.06e16; % Exponential prefactor (Pa^3.5 s)
        A0_wet = 2.05e14; % Exponential prefactor (Pa^3.5 s)
        Ea_dry = 535e3;   % Activation energy (J/mol)
        Ea_wet = 515e3;   % Activation energy (J/mol)
        dTmdP  = 4.85e-8; % Melt temp gradient [K/Pa]
        Tm     = 1478;
        T0     = 1653;
        % Va = (Ea / Tm) * dTm/dP
        Va0    = (Ea_dry/Tm)*dTmdP;
        sigma0 = 0.3e6;   % Reference stess [Pa]
        P0     = 100e6;   % Reference pressure [Pa]
        eta0   = (A0_dry / sigma0^(n-1))*exp(((Ea_dry+Va0*P0)./(MAT.R * T0))); % ref viscosity
        
        y      = 1.15;     % Grüneisen parameters
        
        
        % General
        MAT.SIL.L     = 452e3;  % Latent heat of fusion [kJ/kg]
        MAT.SIL.Tm0   = 1478;   % Reference melting temperature at surface [K]
                                % 2163 K for forsterite
                                % 1240 - 1470 K for basalt (low viscosity)
                                % 930 - 1070 K for Rhyolite (high viscosity)
                                % Tm = 1478 + 4.85e-8 P [T in K, P in Pa]
        
        % Solids
        MAT.SIL.s.nu    = nu;     % Poisson's ratio
        MAT.SIL.s.Emod  = Emod;   % Young's modulus [Pa]
        MAT.SIL.s.Gmod  = Gmod;   % Elastic shear modulus [Pa]
        MAT.SIL.s.Kmod  = Kmod;   % Bulk modulus [Pa]
        MAT.SIL.s.T0    = T0;     % Reference temperature for props [K]
        MAT.SIL.s.eta0  = eta0;   % Reference viscosity at T0 [Pa s]
        MAT.SIL.s.Ea    = Ea_dry; % Activation Energy [J mol-1]
        MAT.SIL.s.Va    = Va0;    % Activation Volume [m^3 mol-1]
        MAT.SIL.s.k0    = 2.47;   % Referemce thermal conductivity [W m-1 K-1]
        MAT.SIL.s.rho0  = 4390;   % Reference density [kg/m^3]
        MAT.SIL.s.Cp0   = 875.5;  % Specific heat capacity [J/kg K]
        
        % Melts
        MAT.SIL.m.Kmod  = 22.6e9; % Bulk modulus [Pa]
        MAT.SIL.m.eta0  = 1e3;    % Viscosity of melt [Pa s]
        MAT.SIL.m.rho0  = 2900;   % Density [kg/m^3]
        MAT.SIL.m.Cp0   = 1600;   % Specific heat capacity [J/kg K]
        
    case 'MB'
        %%%%%
        % Something MORBy (hadgepadge)
        %%%%
        Kmod = 127.9e9;
        Gmod = 50.3e9;
        Emod = 9*Gmod*Kmod/(Gmod+3*Kmod);
        nu   = Emod/(2*Gmod) -1;
        
        n      = 3.5;     % Stress exponent
        d      = 1e-6;    % Grain size [m]
        A0_dry = 2.06e16; % Exponential prefactor (Pa^3.5 s)
        A0_wet = 2.05e14; % Exponential prefactor (Pa^3.5 s)
        Ea_dry = 535e3;   % Activation energy (J/mol)
        Ea_wet = 515e3;   % Activation energy (J/mol)
        dTmdP  = 4.85e-8; % Melt temp gradient [K/Pa]
        Tm     = 1478;
        T0     = 1653;
        % Va = (Ea / Tm) * dTm/dP
        Va0    = (Ea_dry/Tm)*dTmdP;
        sigma0 = 0.3e6;   % Reference stess [Pa]
        P0     = 100e6;   % Reference pressure [Pa]
        eta0   = (A0_dry / sigma0^(n-1))*exp(((Ea_dry+Va0*P0)./(MAT.R * T0))); % ref viscosity
        
        y      = 1.15;     % Grüneisen parameters
        
        
        % General
        MAT.SIL.L     = 452e3;  % Latent heat of fusion [kJ/kg]
        MAT.SIL.Tm0   = 1478;   % Reference melting temperature at surface [K]
                                % 2163 K for forsterite
                                % 1240 - 1470 K for basalt (low viscosity)
                                % 930 - 1070 K for Rhyolite (high viscosity)
                                % Tm = 1478 + 4.85e-8 P [T in K, P in Pa]
        
        % Solids
        MAT.SIL.s.nu    = nu;     % Poisson's ratio
        MAT.SIL.s.Emod  = Emod;   % Young's modulus [Pa]
        MAT.SIL.s.Gmod  = Gmod;   % Elastic shear modulus [Pa]
        MAT.SIL.s.Kmod  = Kmod;   % Bulk modulus [Pa]
        MAT.SIL.s.T0    = T0;     % Reference temperature for props [K]
        MAT.SIL.s.eta0  = eta0;   % Reference viscosity at T0 [Pa s]
        MAT.SIL.s.Ea    = Ea_dry; % Activation Energy [J mol-1]
        MAT.SIL.s.Va    = Va0;    % Activation Volume [m^3 mol-1]
        MAT.SIL.s.k0    = 4.0;    % Referemce thermal conductivity [W m-1 K-1]
        MAT.SIL.s.rho0  = 2900;   % Reference density [kg/m^3]
        MAT.SIL.s.Cp0   = 875.5;  % Specific heat capacity [J/kg K]
        
        % Melts
        MAT.SIL.m.Kmod  = 22.6e9; % Bulk modulus [Pa]
        MAT.SIL.m.eta0  = 1e3;    % Viscosity of melt [Pa s]
        MAT.SIL.m.rho0  = 2400;   % Density [kg/m^3]
        MAT.SIL.m.Cp0   = 1600;   % Specific heat capacity [J/kg K]
        
        
        
    otherwise
        error('Composition for silicates not correctly specified.');
        
end


%%%%%
% Iron
%%%%%
% Find the required iron density to meet the global mass constraint
if isempty(IN.fm0_H2O); IN.fm0_H2O = (4/3)*pi*(BOD.R^3-(BOD.R - IN.Hmax_H2O)^3)*MAT.H2O.s.rho0/BOD.m; end
if isempty(IN.Hmax_H2O); IN.Hmax_H2O = BOD.R-floor((BOD.R^3 - (3/(4*pi))*(IN.fm0_H2O*BOD.m)/MAT.H2O.s.rho0)^(1/3)); end

vIS = BOD.V -  (4/3)*pi*(BOD.R^3-(BOD.R - IN.Hmax_H2O)^3); % Volume of iron + rock

if isempty(IN.fm0_irn)
    IN.fm0_sil  = (4/3)*pi*((BOD.R - IN.Hmax_H2O)^3-IN.Hmax_irn^3)*MAT.SIL.s.rho0/BOD.m;
    IN.fm0_irn  = (1-(IN.fm0_H2O+IN.fm0_sil));
elseif isempty(IN.Hmax_irn)
    IN.fm0_sil  = (1-(IN.fm0_H2O+IN.fm0_irn));
    mSil        = IN.fm0_sil * BOD.m; % Mass of rock
    vSil        = mSil / MAT.SIL.s.rho0; % Volume of rock
    vIrn        = vIS-vSil; % Volume of iron
    IN.Hmax_irn = ((3/(4*pi))*vIrn)^(1/3);
end
rhoIrn = IN.fm0_irn * BOD.m / ((4/3)*pi*IN.Hmax_irn^3);

if rhoIrn<MAT.SIL.s.rho0
    error('Fix core radius/fraction, or silicate material. The core is lighter than the rock');
elseif rhoIrn>9e3
    error('Fix core radius/fraction, or silicate material. The core is over 9000 kg/m^3');
end

% Scratch
Gmod = 34.6e9;
Emod = 2*Gmod*(1+nu);
Kmod = Emod/(3*(1-2*nu));

% General
MAT.IRN.Tm0 = IN.Tm_irn; % Melting temperature [K]
MAT.IRN.L   = 247e3;     % Latent heat of fusion [kJ/kg]

% Solids
MAT.IRN.s.k0    = 27.0;   % Reference Thermal conductivity [W/mK]
MAT.IRN.s.nu    = 0.27;   % Poisson's ratio
MAT.IRN.s.Emod  = Emod;   % Young's modulus [Pa]
MAT.IRN.s.Gmod  = Gmod;   % Elastic shear modulus [Pa]
MAT.IRN.s.Kmod  = Kmod;   % Bulk modulus [Pa]
MAT.IRN.s.eta0  = MAT.IRN.s.Gmod/BOD.omega; % Reference viscosity at T0 [Pa s]
MAT.IRN.s.rho0  = rhoIrn; % Density [kg/m^3]
MAT.IRN.s.Cp0   = 450;    % Specific heat capacity [J/kg K]

% Melts
MAT.IRN.m.Kmod  = 2.2e9;  % Bulk modulus [Pa]
MAT.IRN.m.T0    = 273;    % Reference temperature for props [K]
MAT.IRN.m.eta0  = 1e-3;   % Viscosity of melt  [Pa s]
MAT.IRN.m.rho0  = 0.95*rhoIrn; % Density [kg/m^3]
MAT.IRN.m.Cp0   = 820;   % Specific heat capacity [J/kg K]
















end