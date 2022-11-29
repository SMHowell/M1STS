%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Takes a property on staggered nodes and returns the volumetric average on
% normal nodes
function [prop_n] = s2nMass(M,prop)

rho_s   = n2sVolumetric(M,M.rho); % Staggered density
refIndT = 2:M.Nz -1;              % indices of solution
prop_n  = zeros(1,M.Nz); 
prop_n(refIndT) = (prop(1:end-1).*M.V_s(1:end-1).*rho_s(1:end-1) + prop(2:end).*M.V_s(2:end).*rho_s(2:end))./ ...
                  (M.V_s(1:end-1).*rho_s(1:end-1)+M.V_s(2:end).*rho_s(2:end));

% Set first and last values
prop_n(1)   = prop_n(2);
prop_n(end) = prop_n(end-1);   

end