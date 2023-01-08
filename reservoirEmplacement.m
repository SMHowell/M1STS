%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [COMP, M] = reservoirEmplacement(BOD,IN,M,COMP)

% Check if its time for reservoir emplacement
if (M.t>IN.tRes) && (M.resEmp == 0)

    % Set flag
    M.resEmp = 1;

    % Save temperature field
    Tinit = M.T;
    
    % Set temperature and melt fraction
    M.T(  (M.z  <=IN.zResTop+2*IN.rRes) & (M.z  >=IN.zResTop)) = M.Tm_res;
    M.vfm((M.z_s<=IN.zResTop+2*IN.rRes) & (M.z_s>=IN.zResTop)) = 1;
    
    % Find interfaces
    M.rResTop = M.r(end-1)-IN.zResTop;                % Reservoir top interface radius
    M.iResTop = find((M.rResTop - M.r>0)>0,1,'last'); % Reservoir top interface element index

    M.rResBot = M.r(end-1)-IN.zResTop-2*IN.rRes;      % Reservoir bottom interface radius
    M.iResBot = find((M.rResBot - M.r>0)>0,1,'last'); % Reservoir bottom interface element index
   
    % Set properties
    M.vif    = 0;
    M.vlf    = 1;
    M.vShell = 4/3*pi*M.rResTop^3 - 4/3*pi*M.rResBot^3; % volume of the liquid shell
    M.rRes   = (M.rResTop-M.rResBot)/2; % Reservoir radius
    M.vRes   = (4/3)*pi*M.rRes^3; % Reservoir volume
    M.vRes_old  = M.vRes;
    M.Vice_old  = 0;
    M.vRes_init = M.vRes;         % Initial reservoir volume 
    M.rhoRes = M.rhoOcn;          % Reservoir density
    M.CpRes  = M.CpOcn;           % Reservoir specific heat capacity
    M.mRes   = M.vRes*M.rhoRes;   % Reservoir mass

    % Initialize energy
    [COMP, M] = initializeEnergy(IN,COMP,M);

    % Calculate Maxwell time far from reservoir
    IN.E        = 1e9;                              % Young modulus
    iResCenter  = M.iResBot + round((M.iResTop-M.iResBot)/2);
    M.Tfar      = Tinit(iResCenter);                % far field temperature
    M.eta       = 1e14*exp(25.2*(273/M.Tfar-1));    % viscosity 
    M.tMaxwell  = M.eta / IN.E + M.t;               % Maxwell time (after res. emplacement)   

    % Threshold overpressure before eruption
    sigmaSurf = 1.7e6;    % Ice tensile strength at 100K [Pa]
    sigmaConv = 1e6;      % Ice tensile strength at 250K [Pa]
    sigmaFar  = sigmaSurf+(((sigmaConv-sigmaSurf)/150)*M.Tfar); % far field ice tensile strength
    P0        = IN.zResTop * BOD.g * 920; % Better value for ice density????  <------------ /!\
    M.DeltaPc = 2*(sigmaFar + P0);

end

end
























