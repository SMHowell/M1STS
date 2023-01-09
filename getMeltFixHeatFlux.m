%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = getMeltFixHeatFlux(M,BOD,IN,MAT,qSurfAddtl)

% First, we must get the surface temperature that matches the steady state
% conductive thickness prescribed
sigma = 5.67e-8;     % Stefan Boltzamnn constant [W m^-2 K^-4]
m     = 18.01528e-3; % Molar weight of H2O [kg mol^-1]
R     = 8.314;       % Ideal gas constant [J mol^-1 K^-1]
Es    = 3.8e26;      % Solar output [W]

% Search array
Tsearch = linspace(10,MAT.H2O.Tm0,1000); % Surface temperature range to search [K]

%%%%%%%%%%%%%%%%%%%%%%%
% Emissive heat flux
%%%%%%%%%%%%%%%%%%%%%%%
qEmis = BOD.eps * sigma * Tsearch.^4;


%%%%%%%%%%%%%%%%%%%%%%%
% Solar insolation heat flux
%%%%%%%%%%%%%%%%%%%%%%%
% Oj + Stev https://reader.elsevier.com/reader/sd/pii/0019103589900523
% Average over latitude
switch numel(IN.lat)
    case 1
        coLatEff = IN.lat;
    case 2
        % Find effective latitude. Knowing T~(q_sol)^4 and q_sol~(ob^2+coLat^2)^1/2,
        % we develop a scaling to find the effective average. Away from the
        % equator (where coLat > ob), effective latitude is:
        coLat    = 90 - IN.lat; % Colatitude
        coLatEff = ((4/5) * (coLat(2)^(5/4) - coLat(1)^(5/4)) / (coLat(2)-coLat(1))).^(4); % Effective colatitude
        M.latEff = 90 - coLatEff; % Effective latitude
    otherwise
        % Find effective latitude. Knowing T~(q_sol)^-4 and q_sol~(ob^2+coLat^2)^1/2
        % we develop a scaling to find the effective average. Away from the pole
        % center (where coLat >> ob), effective latitude is:
        coLat    = 90 - [min(IN.lat),max(IN.lat)]; % Colatitude
        coLatEff = ((4/5) * (coLat(2)^(5/4) - coLat(1)^(5/4)) / (coLat(2)-coLat(1))).^(4); % Effective colatitude
        M.latEff = 90 - coLatEff; % Effective latitude
end

Fs0      = Es / (4 * pi * BOD.RParent^2);  % solar insolation at equatiorial high noon [W/m^2]
Fs       = Fs0 * sqrt((BOD.ob*pi/180)^2 + (coLatEff*pi/180).^2) / (pi * sqrt(2)); % Obliquity and colatitude correction
qSol     = (1-BOD.BoloA) * Fs;                 % Average annual solar insolation [W/m^2]


%%%%%%%%%%%%%%%%%%%%%%%
% Sublimation heat flux
%%%%%%%%%%%%%%%%%%%%%%%
% Calculate latent heat of sublimation
% Polynomial curve fits to Table 2.1. R. R. Rogers; M. K. Yau (1989).
% BOD.BoloA Short Course in Cloud Physics (3rd ed.). Pergamon Press. p. 16.
% ISBN 0-7506-3215-1.
L = (2834.1 - 0.29*Tsearch - 0.004*Tsearch.^2)*1e3;

% Calculate the partial pressure of ice
% https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/93GL00105
% NOTE log = log10 in this paper, not ln
pA = -2663.5; % empirical coefficient
pB = 12.537;  % empirical coefficient
P = 10.^(pA./Tsearch + pB); % partial pressure [Pa]

% heat flux
qSub = L .* P .* sqrt(m./(2 * pi * R * Tsearch));


%%%%%%%%%%%%%%%%%%%%%%%
% Required Geologic Heat Flux [W/m^2]
%%%%%%%%%%%%%%%%%%%%%%%
qIn   = qSol + qSurfAddtl; % Paid into heat budget
qOut  = qEmis + qSub; % Withdrawn from heat budget
qAnom = qOut-qIn;     % Anomalous heat flux



%%%%%%%%%%%%%%%%%%%%%%%
% Find Nominal Surface Temp
%%%%%%%%%%%%%%%%%%%%%%%
Tsurf_0 = M.Tsurf; 
qIce    = qSurfAddtl;

z     = linspace(0,IN.Hmin_H2O,1000); % Depth profile [m]
H_ice = IN.Hmin_H2O; % Layer thickness

% Set lower temperature
Tb = MAT.H2O.Tm0;

dT_crit = 1e-5; % Convergence criterion for surface temperature [K]
T_old   = 0;    % Old surface temperature [K]
while abs(Tsurf_0-T_old) > dT_crit
    T_z   = Tsurf_0*(Tb/Tsurf_0).^(z/H_ice);             % Temperature profile
    k     = mean(632 ./ T_z + 0.38 - 1.97e-3 .* T_z);    % Mean solid ice conductivity [W/m K]
    qIce  = (k * (Tb - Tsurf_0)) * (1/H_ice + 1/BOD.R);  % Heat flux through ice evaluated at the surface [W/m^2]
    T_old = Tsurf_0;                                     % Store previous result
    Tsurf_0  = interp1(qAnom-qIce,Tsearch,0);            % New surface temperature guess
end  

%%%%%%%%%%%%%%%%%%%%%%%
% Store Results
%%%%%%%%%%%%%%%%%%%%%%%
M.Tsurf   = Tsurf_0; % Surface temperature [K]
M.T(end)  = M.Tsurf;

M.qSol    = qSol; % Solar insolation heat flux
M.qLoss   = qIce; % Heat flux from interior through ice shell [W/m^2]

% Enforce space BCs
M.T(end-1) = M.Tsurf;


end








