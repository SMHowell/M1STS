%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% reservoir freezing and evolution.
% Sam Howell, Elodie Lesage
% samuel.m.howell@jpl.nasa.gov
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = reservoirEmplacement(M,IN)

% Check if its time for reservoir emplacement
if (M.t>IN.tRes) && (M.resEmp == 0)
    % Set flag
    M.resEmp = 1;
    
    % Set temperature and melt fraction
    M.T(  (M.z  <=IN.zResTop+2*IN.rRes) & (M.z  >=IN.zResTop)) = M.Tm_res;
    M.vfm((M.z_s<=IN.zResTop+2*IN.rRes) & (M.z_s>=IN.zResTop)) = 1;
    
    % Find interfaces
    M.rResTop = M.r(end-1)-IN.zResTop;                % Reservoir top interface radius
    M.iResTop = find((M.rResTop - M.r>0)>0,1,'last'); % Reservoir top interface element index

    M.rResBot = M.r(end-1)-IN.zResTop-2*IN.rRes;      % Reservoir bottom interface radius
    M.iResBot = find((M.rResBot - M.r>0)>0,1,'last'); % Reservoir bottom interface element index
   
    % Set properties
    M.rRes   = (M.rResTop-M.rResBot)/2; % Reservoir radius
    M.vRes   = (4/3)*pi*M.rRes^3; % Reservoir volume
    M.rhoRes = M.rhoOcn;          % Reservoir density
    M.CpRes  = M.CpOcn;           % Reservoir specific heat capacity
    M.mRes   = M.vRes*M.rhoRes;   % Reservoir mass

end

end
























