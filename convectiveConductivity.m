%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k_A, k_B] = convectiveConductivity(M,k_A,k_B)
% Set thermal conductivity on interface elements for the ocean
vfs_intf = (1-M.vfm(M.iOcnTop-1:M.iOcnTop+1)); % Solid volume fraction of elements
V_intf   = M.V_s(M.iOcnTop-1:M.iOcnTop+1);     % Element volumes
k_intf   = (sum(V_intf.*vfs_intf./M.k(M.iOcnTop-1:M.iOcnTop+1))./ sum(V_intf))^(-1); % Effective conductivity

k_A(M.iOcnTop-1:M.iOcnTop+1) = k_intf;
k_B(M.iOcnTop-2:M.iOcnTop)   = k_intf;

% Set thermal conductivity on interface elements for the reservoir top and
% bottom if they have been prescribed
if M.vRes>0      
    % Reservoir top
    vfs_intf = (1-M.vfm(M.iResTop-1:M.iResTop+1)); % Solid volume fraction of elements
    V_intf   = M.V_s(M.iResTop-1:M.iResTop+1);     % Element volumes
    k_intf   = (sum(V_intf.*vfs_intf./M.k(M.iResTop-1:M.iResTop+1))./ sum(V_intf))^(-1); % Effective conductivity

    k_A(M.iResTop-1:M.iResTop+1) = k_intf;
    k_B(M.iResTop-2:M.iResTop)   = k_intf;
    
    % Reservoir bottom
    vfs_intf = (1-M.vfm(M.iResBot-1:M.iResBot+1)); % Solid volume fraction of elements
    V_intf   = M.V_s(M.iResBot-1:M.iResBot+1);     % Element volumes
    k_intf   = (sum(V_intf.*vfs_intf./M.k(M.iResBot-1:M.iResBot+1))./ sum(V_intf))^(-1); % Effective conductivity

    k_A(M.iResBot-1:M.iResBot+1) = k_intf;
    k_B(M.iResBot-2:M.iResBot)   = k_intf;
end


end



















