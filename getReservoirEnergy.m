%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = getReservoirEnergy(M,Q_A,Q_B)

%%%%%%%%%%%%%%%%%%%%%%%
% Initialize thermal properties
%%%%%%%%%%%%%%%%%%%%%%%
% Track net heat to ocean, remembering that the first entry in Q_A and Q_B
% corresponds to the second entry in T. Thus, we will subtract 1 from each
% of the indices, and find the heat leaving the seafloor and entering the
% ice shell.
dEocn_top = 4*pi*M.dt*(Q_B(M.iOcnTop-1)); % Energy out of ocean from ice
dEocn_bot = 4*pi*M.dt*(Q_A(M.iOcnBot));   % Energy into ocean from rock
M.dE_ocn  = dEocn_bot-dEocn_top;          % Net energy change in ocean

M.dEres_bot = 4*pi*M.dt*(Q_A(M.iResBot));    % Energy into reservoir from bottom
M.dEres_top = -4*pi*M.dt*(Q_B(M.iResTop-1)); % Energy out of reservoir from top

end










































