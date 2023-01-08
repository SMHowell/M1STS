%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T, rho, vlf, vif] = coolFreeze(M, COMP, IN)

% find index of closest remaining energy fraction from composition data
i_low = find(COMP.fE_rmn{IN.simu} >= M.fE_rmn, 1, 'last');
i_max = length(COMP.fE_rmn{IN.simu});

% if i_low == 3
%     disp(i)
% end

% interpolate values of temperature, volumic frozen fraction and density
% associated
if i_low == 1 
    indexInterp = [i_low, i_low+1];
elseif i_low+1 == i_max
    indexInterp = [i_low-1, i_low, i_low+1];
elseif i_low == i_max
    indexInterp = [i_low-2, i_low-1, i_low];
else
    indexInterp = [i_low-1, i_low, i_low+1, i_low+2];
end

% interpolation:
if i_low == 1   % linear interpolation between 1st and second data points
    T           = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.temps{IN.simu}(indexInterp),M.fE_rmn,'linear');
    rho         = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.rho_l{IN.simu}(indexInterp),M.fE_rmn,'linear');
    vlf         = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.vlf{IN.simu}(indexInterp),M.fE_rmn,'linear');
    vif         = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.vif{IN.simu}(indexInterp),M.fE_rmn,'linear');
else
    T           = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.temps{IN.simu}(indexInterp),M.fE_rmn,'spline');
    rho         = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.rho_l{IN.simu}(indexInterp),M.fE_rmn,'spline');
    vlf         = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.vlf{IN.simu}(indexInterp),M.fE_rmn,'spline');
    vif         = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.vif{IN.simu}(indexInterp),M.fE_rmn,'spline');
end

end