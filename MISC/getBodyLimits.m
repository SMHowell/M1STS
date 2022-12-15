%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric thermal solver
% Conservative finite differences, explicit forward in time
% Samuel.M.Howell@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%

function [LIMS] = getBodyLimits(BC,rCore)

% Prep
eutectic = load('EutecticData/FeS.mat'); % Core melting properties

% Body parameters (Ganymede)
m_body  = BC.m_body; % Body mass [kg]
r_body  = BC.r_body; % Body radius [km]
r_sil   = BC.r_body - BC.D_H2O; % Mantle outer radius [km]
rho_H2O = BC.rho_H2O; % H2O layer density  [kg/m^3]

% Universal params
G = 6.6743e-11; % UGC

% Calculate layer values for easy book keeping
r_core = rCore;                   % Radius used in integration [km]
V_body = (4/3)*pi*r_body^3;       % Volume body [m^3]
V_bulk = (4/3)*pi*r_sil^3;        % Volume of non-H2O stuff [m^3]
V_H2O  = V_body - V_bulk;         % Volume of H2O layer [m^3]
V_core = (4/3)*pi*r_core.^3;      % Volume of core [m^3]
V_sil  = V_bulk - V_core;         % Volume of silicates [m^3] 

m_H2O  = V_H2O * rho_H2O;    % Mass of water [kg]
m_bulk = m_body - m_H2O;     % Mass of rock and metal [kg]

% Density physical limits
rho_sil_lim_min  = 2500;     % Mantle density lower bound
rho_core_lim_min = 5000; % Lower core desnsity limit [kg/m^3]
rho_core_lim_max = 8000; % Upper core desnsity limit [kg/m^3]

% Set core density bounds
LIMS.rho_core_min = repmat(rho_core_lim_min,1,size(r_core,2));
LIMS.rho_core_max = repmat(rho_core_lim_max,1,size(r_core,2));

% Check for core densities that violate silicate limit
LIMS.rho_core_max = min(LIMS.rho_core_max,(m_body - m_H2O - rho_sil_lim_min.*V_sil)./V_core);

% Calculate mantle density bounds [kg/m^3]
LIMS.rho_sil_min = (m_body - m_H2O - V_core.*LIMS.rho_core_max)./(V_sil);
LIMS.rho_sil_max = (m_body - m_H2O - V_core.*LIMS.rho_core_min)./(V_sil);


%% For each body, calculate the pressure as a function of depth and get temperature ranges
r_int  = linspace(0,r_body,1000); % Radius used in integration [km]
rho_hi = zeros(size(r_int)); % Initialize body density upper bound [kg/m^3]
rho_lo = zeros(size(r_int)); % Initialize body density lower bound [kg/m^3]

P_cmb_hi = zeros(size(r_core)); % Lower pressure bounds
P_cmb_lo = zeros(size(r_core)); % Upper pressure bounds
LIMS.T_melt_lo = zeros(size(r_core)); % Lower temp bounds
LIMS.T_melt_hi = zeros(size(r_core)); % Upper temp bounds
for i = 1:length(r_core)
    r_core_int = r_core(i);
    
    % Get mantle density bounds
    rho_m_lo = interp1(r_core,LIMS.rho_sil_min,r_core_int);
    rho_m_hi = interp1(r_core,LIMS.rho_sil_max,r_core_int);

    % Get mantle core bounds
    rho_c_lo = interp1(r_core,LIMS.rho_core_min,r_core_int);
    rho_c_hi = interp1(r_core,LIMS.rho_core_max,r_core_int);
    
    % Get denisty lower bound
    rho_lo(r_int<r_core_int) = rho_c_lo;
    rho_lo(r_int>=r_core_int & r_int<=r_sil) = rho_m_lo;
    rho_lo(r_int>r_sil) = rho_H2O;
    
    % Get density upper bounds
    rho_hi(r_int<r_core_int) = rho_c_hi;
    rho_hi(r_int>=r_core_int & r_int<=r_sil) = rho_m_hi;
    rho_hi(r_int>r_sil) = rho_H2O;
    
    % Get mass derivative
    dmdr_lo = 4*pi*r_int.^2.*rho_lo;
    dmdr_hi = 4*pi*r_int.^2.*rho_hi;
    
    % Get mass
    m_lo = cumtrapz(r_int,dmdr_lo);
    m_hi = cumtrapz(r_int,dmdr_hi);
    
    % Get pressure derivative
    dPdr_lo = -G * m_lo .* rho_lo ./ (r_int.^2+1e-5);
    dPdr_hi = -G * m_hi .* rho_hi ./ (r_int.^2+1e-5);
    
    % Integrate for pressure
    P_lo = -cumtrapz(r_int,dPdr_lo);
    P_hi = -cumtrapz(r_int,dPdr_hi);
    
    % Corection to match Schubert, 2009
    P_lo = max(P_lo)-P_lo;
    P_hi = max(P_hi)-P_hi;
    
    % Find core pressure bounds
    P_cmb_lo(i) = interp1(r_int,P_lo,r_core_int);
    P_cmb_hi(i) = interp1(r_int,P_hi,r_core_int);
    
    % Now calculate eutectic bounds
    FeS_1 = interp1(eutectic.FeS_data{1}(:,1)*1e6,eutectic.FeS_data{1}(:,2),P_cmb_hi(i));
    if isnan(FeS_1)
        FeS_1 = max(eutectic.FeS_data{1}(:,2));
    end
    
    FeS_2 = interp1(eutectic.FeS_data{2}(:,1)*1e6,eutectic.FeS_data{2}(:,2),P_cmb_hi(i));
    FeS_3 = interp1(eutectic.FeS_data{3}(:,1)*1e6,eutectic.FeS_data{3}(:,2),P_cmb_hi(i));
    FeS_4 = interp1(eutectic.FeS_data{4}(:,1)*1e6,eutectic.FeS_data{4}(:,2),P_cmb_hi(i));
    
    FeS_5 = interp1(eutectic.FeS_data{1}(:,1)*1e6,eutectic.FeS_data{1}(:,2),P_cmb_lo(i));
    if isnan(FeS_1)
        FeS_1 = max(eutectic.FeS_data{1}(:,2));
    end
    FeS_6 = interp1(eutectic.FeS_data{2}(:,1)*1e6,eutectic.FeS_data{2}(:,2),P_cmb_lo(i));
    FeS_7 = interp1(eutectic.FeS_data{3}(:,1)*1e6,eutectic.FeS_data{3}(:,2),P_cmb_lo(i));
    FeS_8 = interp1(eutectic.FeS_data{4}(:,1)*1e6,eutectic.FeS_data{4}(:,2),P_cmb_lo(i));
    
    % Save bounds
    LIMS.T_melt_lo(i) = nanmin([FeS_1,FeS_2,FeS_3,FeS_4,FeS_5,FeS_6,FeS_7,FeS_8]);
    LIMS.T_melt_hi(i) = nanmax([FeS_1,FeS_2,FeS_3,FeS_4,FeS_5,FeS_6,FeS_7,FeS_8]);
    
end


end















