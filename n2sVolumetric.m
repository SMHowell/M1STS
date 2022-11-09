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
function [prop_s] = n2sVolumetric(M,prop)

f1 = (M.r_s(M.refInd).^3 - M.r(M.refInd).^3)./  (M.r(M.refInd+1).^3 - M.r(M.refInd).^3); % Weighting for element i
f2 = (M.r(M.refInd+1).^3 - M.r_s(M.refInd).^3)./(M.r(M.refInd+1).^3 - M.r(M.refInd).^3); % Weighting for element i+1

prop_s = f1 .* prop(M.refInd) + f2 .* prop(M.refInd+1);

end