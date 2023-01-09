%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k_A, k_B] = boundaryConductivity(M,k_A,k_B)
% Set thermal conductivity on interface elements for the ocean
if M.vOcn > 0 && M.iOcnTop < M.Nz-2
    vfm = M.mat.fm_s(M.mat.iH2Omelt,:);
        
    vfs_intf = (1-vfm(M.iOcnTop-1:M.iOcnTop+1)); % Solid volume fraction of elements
    V_intf   = M.V_s(M.iOcnTop-1:M.iOcnTop+1);   % Element volumes    
    k_intf   = (sum(V_intf.*vfs_intf./M.k(M.iOcnTop-1:M.iOcnTop+1))./ sum(V_intf))^(-1); % Effective conductivity
    
    k_A(M.iOcnTop-1:M.iOcnTop+1) = k_intf;
    k_B(M.iOcnTop-2:M.iOcnTop)   = k_intf;
    
    k_A(M.iOcnBot-1:M.iOcnBot+1) = k_intf;
    k_B(M.iOcnBot-2:M.iOcnBot)   = k_intf;
end


end



















