%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric thermal solver
% Conservative finite differences, explicit forward in time
% Samuel.M.Howell@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = getIceRheology(M,BOD,MAT)

% Inputs for aggregate considerations
vfmCr      = 0.40;   % Melt fraction required for disaggregation (~0.4 - 0.6):
                     % https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2002JE001943
aMelt      = 30;     % Exponent on melt-fraction dependence of viscosity: 
                     % soest.hawaii.edu/GG/FACULTY/smithkonter/GG631/other/HirthKohlstedt_2000.pdf

% Material Inputs 
MAT.H2O.s.eta0 = MAT.H2O.s.Gmod/BOD.omega; % Reference viscosity at MAT.H2O.s.T0 [Pa s]

% Get creep viscosity
etaDiff = MAT.H2O.s.eta0 * exp((MAT.H2O.s.Ea/(MAT.R*MAT.H2O.s.T0))*(MAT.H2O.s.T0./M.T-1));

% Get brittle viscosity limit
Z       = MAT.H2O.s.Gmod*M.dt./(MAT.H2O.s.Gmod*M.dt+etaDiff); % Elasticity correction
M.etaVE = Z.*etaDiff;

% Consider the role of *partial* melt!
iceInd   = (M.iOcnTop+2:M.Nz-1);    % Indices contained entirely w/in the ice shell
fVmelt_n = s2nVolumetric(M,M.mat.fV_s(M.mat.iH2Omelt,:));
etaIce   = M.etaVE(iceInd).*exp(-aMelt*fVmelt_n(iceInd));

% % Einstein Roscoe relationship above BOD.vfmCr
% etaIce(fVmelt_n>=vfmCr) = MAT.H2O.m.eta0.*(1.35*fVmelt_n(fVmelt_n>vfmCr)-0.35).^(-5/2);

% Composite
M.etaVE(iceInd) = etaIce;


end