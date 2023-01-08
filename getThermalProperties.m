%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = getThermalProperties(M)

% Get thermal properties assuming solid ice for nodes
% Constants
T0    = 270;      % Ref temp for laws [K]
a0    = 1.704e-4; % Ref thermal expansivity [1/K]
rho0  = 917;      % Ref density @ T0 [kg/m^3]

% Density
M.alpha = a0 * M.T / T0;                           % Thermal expansivity [1/K]
M.rho   = (1-M.phi).*rho0.*(1-(M.T-T0).*M.alpha);  % Density [kg/m^3]

% Conductivity
k0_c = 2.107; % Ref thermal cond of continuous phase @ T0 [W/ m K]
k0_d = 0;     % Ref thermal cond of dispersed phase @ T0 [W/ m K]
kT_c = k0_c*T0./M.T;   % Continuous phase conductivity [W/m K]
kT_d = k0_d;   % Dispersed phase conductivity [W/m K]

k_s_fp = kT_c .* (M.phi - 1).^2 + kT_d; % Snow
k_f_fp = kT_c .* ((1-M.phi) + M.phi .* (3 * k0_d) ./ (2.*kT_c + kT_d))./((1-M.phi) + M.phi .* (3 * kT_c) ./ (2*kT_c + kT_d)); % Firn
M.k    = (M.phi./k_s_fp + (1-M.phi)./(k_f_fp)).^-1; % Composite

% Heat capacity
M.Cp =  1e3 * (7.73e-3 * M.T .* (1 - exp(-1.263e-3 * M.T.^2)) .* ... % Specific heat capacity [J/kg K] 10.1051/0004-6361:20031746
           (1 + exp(-3 * M.T.^(1/2)) * 8.47e-3 .* M.T.^6 + ...
            2.0825e-7 * M.T.^4 .* exp(-4.97e-2 * M.T))); 

%%%%%%%%
% Account for melt contribution to Cp, rho, but *NOT* k, which is handled in
% the thermal solver and convectiveConductivity.m. 
% Lay out ocean and reservoir thermal properties.
%%%%%%%%

rhoMelt = zeros(size(M.rho));
rhoMelt(M.r<=M.rOcnTop) = M.rhoOcn;
if M.vRes>0; rhoMelt((M.r>= M.rResBot) & (M.r<=M.rResTop)) = M.rhoRes; end

CpMelt = zeros(size(M.Cp));
CpMelt(M.r<=M.rOcnTop) = M.CpOcn;
if M.vRes>0; CpMelt(M.r>= M.rResBot & M.r<=M.rResTop) = M.CpRes; end


% Get Density via volume melt fraction
vfm     = s2nVolumetric(M,M.vfm); % Melt fraction on nodes
M.rho   = M.rho .* (1-vfm) + rhoMelt .* vfm; % Density

% Get Cp via mass melt fraction
Nz      = numel(M.r);
rho_s   = n2sVolumetric(M,M.rho); % Staggered density
refIndT = 2:Nz-1;             % Indices of solution
mfm     = zeros(1,Nz);         % Mass melt fraction
mfm(refIndT) = (M.vfm(1:end-1).*M.V_s(1:end-1).*rho_s(1:end-1) + ... 
                M.vfm(2:end).*M.V_s(2:end).*rho_s(2:end))./ ...
               (M.V_s(1:end-1).*rho_s(1:end-1)+M.V_s(2:end).*rho_s(1:end-1));
M.Cp  = M.Cp .* (1-mfm) + CpMelt .* mfm; % Soecific heat capacity
        
% Diffusivity
M.K = M.k./(M.rho .* M.Cp);


end



















