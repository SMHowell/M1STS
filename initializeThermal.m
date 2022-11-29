%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IN, M] = initializeThermal(IN,BOD,M)

% First, we must get the surface temperature that matches the steady state
% conductive thickness prescribed
sigma  = 5.67e-8;     % Stefan Boltzamnn constant [W m^-2 K^-4]
m      = 18.01528e-3; % Molar weight of H2O [kg mol^-1]
R      = 8.314;       % Ideal gas constant [J mol^-1 K^-1]
Es     = 3.8e26;      % Solar output [W]

% Search array
TsurfTemp = linspace(10,273,1000); % Surface temperature range to search [K]


%%%%%%%%%%%%%%%%%%%%%%%
% Radiogenic heat flux
%%%%%%%%%%%%%%%%%%%%%%%
qRad = BOD.m * 4.5e-12 / (4 * pi * BOD.R^2); % Radiogenic heat flux [W/m^2]


%%%%%%%%%%%%%%%%%%%%%%%
% Emissive heat flux
%%%%%%%%%%%%%%%%%%%%%%%
qEmis = BOD.eps * sigma * TsurfTemp.^4;


%%%%%%%%%%%%%%%%%%%%%%%
% Solar insolation heat flux
%%%%%%%%%%%%%%%%%%%%%%%
% Oj + Stev https://reader.elsevier.com/reader/sd/pii/0019103589900523
Fs0  = Es / (4 * pi * BOD.RParent^2);   % solar insolation at equatiorial high noon [W/m^2]
clat = 90 - IN.lat;                         % Co-latitude
Fs   = Fs0 * sqrt((BOD.ob*pi/180)^2 + (clat*pi/180)^2) / (pi * sqrt(2)); % Obliquity and colatitude correction
qSol = (1-BOD.A) * Fs; % Average annual solar insolation [W/m^2]


%%%%%%%%%%%%%%%%%%%%%%%
% Sublimation heat flux
%%%%%%%%%%%%%%%%%%%%%%%
% Calculate latent heat of sublimation
% Polynomial curve fits to Table 2.1. R. R. Rogers; M. K. Yau (1989).
% BOD.A Short Course in Cloud Physics (3rd ed.). Pergamon Press. p. 16.
% ISBN 0-7506-3215-1.
L = (2834.1 - 0.29*TsurfTemp - 0.004*TsurfTemp.^2)*1e3;

% Calculate the partial pressure of ice
% https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/93GL00105
% NOTE log = log10 in this paper, not ln
pA = -2663.5; % empirical coefficient
pB = 12.537;  % empirical coefficient
P = 10.^(pA./TsurfTemp + pB); % partial pressure [Pa]

% heat flux
qSub = L .* P .* sqrt(m./(2 * pi * R * TsurfTemp));


%%%%%%%%%%%%%%%%%%%%%%%
% Required Geologic Heat Flux [W/m^2]
%%%%%%%%%%%%%%%%%%%%%%%
qIn   = qSol;          % Paid into heat budget
qOut  = qEmis + qSub; % Withdrawn from heat budget
qAnom = qOut-qIn;     % Anomalous heat flux


%%%%%%%%%%%%%%%%%%%%%%%
% Find Nominal Surface Temp
%%%%%%%%%%%%%%%%%%%%%%%
Tsurf_0_init  = interp1(qAnom,TsurfTemp,0); % First guess without geologic heatflow [K]
Tsurf_0       = Tsurf_0_init;              % Initial Surface temperature [K]
z             = linspace(0,IN.H0_ice,1000); % Depth profile [m]

dT_crit = 1e-5; % Convergence criterion for surface temperature [K]
T_old   = 0;    % Old surface temperature [K]
while abs(Tsurf_0-T_old) > dT_crit
    T_z   = Tsurf_0*(IN.Tm_ocn/Tsurf_0).^(z/IN.H0_ice);            % Temperature profile
    k     = mean(632 ./ T_z + 0.38 - 1.97e-3 .* T_z);          % Mean solid ice conductivity [W/m K]
    qIce  = (k * (IN.Tm_ocn - Tsurf_0)) * (1/IN.H0_ice + 1/BOD.R); % Heat flux through ice evaluated at the surface [W/m^2]
    T_old = Tsurf_0;                                           % Store previous result
    Tsurf_0  = interp1(qAnom-qIce,TsurfTemp,0);                % New surface temperature guess
end


%%%%%%%%%%%%%%%%%%%%%%%
% Store Results
%%%%%%%%%%%%%%%%%%%%%%%
M.Tsurf_0 = Tsurf_0;   % Initial surface temperature [K]
M.Tsurf   = M.Tsurf_0; % Surface temperature [K]

M.QIce_0  = qIce*4*pi*BOD.R^2; % Initial required heat flow to balance T (will be applied to base of shell) [W]
M.qSol    = qSol; % Solar insolation heat flux
M.qRad    = qRad; % Radiogenic heat flux
M.qLoss   = qIce; % Heat flux from interior through ice shell [W/m^2]


%%%%%%%%%%%%%%%%%%%%%%%
% Now, thermal structure can be created
%%%%%%%%%%%%%%%%%%%%%%%
M.T = zeros(1,M.Nz);

% Enforce space BCs
M.T(end) = M.Tsurf_0;


%%%%%%%%%%%%%%%%%%%%%%%
% ICE
%%%%%%%%%%%%%%%%%%%%%%%
% Steady state conductive structure
kT = 1;         % 0 -> constant k, 1 -> k ~ 1/T
r1 = M.rOcn;    % Inner radius
r2 = BOD.R;     % Outer radius

T1 = IN.Tm_ocn; % Inner temeprature
T2 = M.Tsurf;   % Outer temperature

t_ind  = find(M.r>=r1 & M.r<=r2); % Temporary indexing
r_temp = M.r(t_ind);
gamma  = ((r2./r_temp).*(r_temp-r1)/(r2-r1)); % Temperature scaling

% Spherically corrected conductive profile
if kT
    % k ~1/T
    T = T1*(T2/T1).^gamma; % Spherically corrected conductive profile
else
    % k is constant
    T = T1+(T2-T1).*gamma;
end

M.T(t_ind) = T;


%%%%%%%%%%%%%%%%%%%%%%%
% Ocean
%%%%%%%%%%%%%%%%%%%%%%%
r1 = M.rSil; % Inner radius
r2 = M.rOcn;  % Outer radius
t_ind  = find(M.r>=r1 & M.r<=r2); % Temporary indexing

% Enforce ocean
M.T(t_ind) = IN.Tm_ocn;


%%%%%%%%%%%%%%%%%%%%%%%
% ROCK
%%%%%%%%%%%%%%%%%%%%%%%
% Steady state conductive structure
kT = 0;         % 0 -> constant k, 1 -> k ~ 1/T
r1 = M.rIrn;    % Inner radius
r2 = M.rSil;    % Outer radius

T1 = IN.T0_irn; % Inner temeprature
T2 = IN.Tm_ocn; % Outer temperature

t_ind  = find(M.r>=r1 & M.r<=r2); % Temporary indexing
r_temp = M.r(t_ind);
gamma  = ((r2./r_temp).*(r_temp-r1)/(r2-r1)); % Temperature scaling

% Spherically corrected conductive profile
if kT
    % k ~1/T
    T = T1*(T2/T1).^gamma; % Spherically corrected conductive profile
else
    % k is constant
    T = T1+(T2-T1).*gamma;
end

M.T(t_ind) = T;


%%%%%%%%%%%%%%%%%%%%%%%
% IRON
%%%%%%%%%%%%%%%%%%%%%%%
% Steady state conductive structure
kT = 0;         % 0 -> constant k, 1 -> k ~ 1/T
r1 = 0;         % Inner radius
r2 = M.rIrn;    % Outer radius

T1 = IN.T0_irn; % Inner temeprature
T2 = IN.T0_irn; % Outer temperature

t_ind  = find(M.r>=r1 & M.r<=r2); % Temporary indexing
r_temp = M.r(t_ind);
gamma  = ((r2./r_temp).*(r_temp-r1)/(r2-r1)); % Temperature scaling

% Spherically corrected conductive profile
if kT
    % k ~1/T
    T = T1*(T2/T1).^gamma; % Spherically corrected conductive profile
else
    % k is constant
    T = T1+(T2-T1).*gamma;
end

M.T(t_ind) = T;
M.T(1)     = M.T(2); % dT/dr -> 0 at r = 0


%%%%%%%%%%%%%%%%%%%%%%%
% Initialize melts
%%%%%%%%%%%%%%%%%%%%%%%
% *volume* fraction (on elements)
M.vfm   = zeros(1,M.Nz-1);         
M.vfm(M.r>=M.rSil & M.r< M.rOcn) = 1; % Impose ocean
M.dm_dt = 0;  % Change in melt mass vs time


%%%%%%%%%%%%%%%%%%%%%%%
% Initialize reservoir/ocean parameters
%%%%%%%%%%%%%%%%%%%%%%%
M.fV     = 0; % Frozen volume fraction of reservoir
M.resEmp = 0; % Flag for empalcement of reservoir

% Track element containing interface
M.iOcnTop = find((M.vfm>0 & M.r_s >= M.rSil & M.r_s <= M.rOcn),1,'last'); % Ocean top interface element index
M.iOcnBot = find((M.vfm>0 & M.r_s >= M.rSil & M.r_s <= M.rOcn),1,'first'); % Ocean bottom interface element index

M.iResTop = []; % Reservoir top interface element index
M.rResTop = []; % Reservoir top interface radius

M.iResBot = []; % Reservoir bottom interface element index
M.rResBot = []; % Reservoir bottom interface radius

M.rRes = []; % Reservoir radius
M.vRes = []; % Reservoir volume
M.mRes = []; % Reservoir mass


%%%%%%%%%%%%%%%%%%%%%%%
% Get Thermal Properties
%%%%%%%%%%%%%%%%%%%%%%%
M.k   = zeros(1,M.Nz);
M.rho = zeros(1,M.Nz);
M.Cp  = zeros(1,M.Nz);

[M] = getPressure(M,BOD);
[M] = getThermalProperties(M,IN);


%%%%%%%%%%%%%%%%%%%%%%%
% Nusselt Number on Elements
%%%%%%%%%%%%%%%%%%%%%%%
M.Ra_cr = 0; % Critical Rayleigh Number
M.Nu    = ones(1,M.Nz-1);


end
























