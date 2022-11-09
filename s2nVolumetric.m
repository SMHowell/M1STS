%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% reservoir freezing and evolution.
% Sam Howell, Elodie Lesage
% samuel.m.howell@jpl.nasa.gov
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Takes a property on staggered nodes and returns the volumetric average on
% normal nodes
function [prop_n] = s2nVolumetric(M,prop)

Nz = numel(M.r);
refIndT = 2:Nz -1; % indices of solution
prop_n = zeros(1,Nz); 
prop_n(refIndT) = (prop(1:end-1).*M.V_s(1:end-1) + prop(2:end).*M.V_s(2:end))./ ...
                  (M.V_s(1:end-1)+M.V_s(2:end));

end