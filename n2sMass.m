%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Takes a property on normal nodes and returns the volumetric average on
% staggered elements
function [prop_s] = n2sMass(M,prop)

rho_n   = M.rho;     % Node density
refIndT = 2:M.Nz-1;  % indices of solution
prop_s  = zeros(1,M.Nz-1); 
prop_s(refIndT) = (prop(1:end-1).*M.V(1:end-1).*rho_n(1:end-1) + prop(2:end).*M.V(2:end).*rho_n(2:end))./ ...
                  (M.V(1:end-1).*rho_n(1:end-1)+M.V(2:end).*rho_n(2:end));

% Set first and last values
prop_s(1)   = prop(1);
prop_s(end) = prop(end);   

end