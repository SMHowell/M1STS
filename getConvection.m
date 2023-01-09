%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = getConvection(M,MAT)

% **** NOTE: REPLACE *****

%%%%%%%%%%%%%%%%%%%%%%%
% Rayleigh Number
%%%%%%%%%%%%%%%%%%%%%%%
M.Ra_cr = 0; 

%%%%%%%%%%%%%%%%%%%%%%%
% Nusselt Number on Elements
%%%%%%%%%%%%%%%%%%%%%%%
M.Nu = ones(1,M.Nz-1);

if M.t> 1.5*3.15e16 % x Gyr
    fV_n = s2nVolumetric(M,M.mat.fV_s);
    M.Nu = 5 * fV_n(M.mat.iSilSolid,:);
end

end
























