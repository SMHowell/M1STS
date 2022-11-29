%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = meltingFreezing(M,IN)

%%%%%%%%%%%%%%%%%%
% Ocean Melting/Freezing
%%%%%%%%%%%%%%%%%%

% Check whether the interface has warmed above Tm
dEice_top = (M.T(M.iOcnTop+1)-IN.Tm_ocn)*M.rho(M.iOcnTop+1)*M.Cp(M.iOcnTop+1)*(4/3)*pi*(M.r_s(M.iOcnTop+1)^3 - M.rOcn^3);
dEice_top(dEice_top<0) = 0;         % Ice is allowed to be colder than melting
if M.T(M.iOcnTop+1)>IN.Tm_ocn % Ice is not allowed to be warmer, energy here is tracked already
    M.T(M.iOcnTop+1) = IN.Tm_ocn;
end

% Calculate total heat gained by ocean
dE = M.dE_ocn + dEice_top;

% Calculate change in melt
dm = dE/M.LH2O;   % Total change in melt mass
dv = dm/M.rhoOcn; % Total change in melt volume

M.dm_dt = dm/M.dt;

% New interface location
M.rOcn = ((3/(4*pi))*dv+M.rOcn^3)^(1/3);
M.iOcnTop = find((M.rOcn - M.r>0)>0,1,'last');

% Locate interfaces
M.iOcnTop = find((M.rOcn - M.r>0)>0,1,'last'); % Ocean top interface element index
M.iOcnBot = find((M.rSil - M.r>0)>0,1,'last'); % Ocean bottom interface element index

% Distribute melt fraction (needed for thermal properties)
if (M.rOcn <= M.rSil) && (abs(dE)>0)
    % Fully frozen ocean
    % Energy Accounting
    rho_s = n2sVolumetric(M,M.rho);
    Cp_s  = n2sVolumetric(M,M.Cp);

    dE = (4/3)*pi*(M.rOcn^3 - M.rSil^3)*M.rhoOcn*M.LH2O;
    dT = dE/(M.V_s(M.iOcnTop)*rho_s(M.iOcnTop)*Cp_s(M.iOcnTop));
    M.T(M.iOcnTop) = M.T(M.iOcnTop)+dT;

    % Close out frozen ocean
    if isempty(M.rResBot) % No melts in the ice
        M.vfm(M.r>=M.rSil) = 0;
    else
        M.vfm(M.r>=M.rSil & M < M.rResBot) = 0;
    end
    M.iOcnBot = M.iOcnBot;
    M.rOcn    = M.rSil;
else
    % Nominal case
    % Set melt fractions
    M.vfm(M.r_s>M.rSil & M.r_s<M.rOcn) = 1;
    M.vfm(M.r_s>M.rOcn | M.r_s<M.rSil) = 0;
    M.vfm(M.iOcnTop) = (4/3)*pi*(M.rOcn^3 - M.r(M.iOcnTop)^3)/M.V_s(M.iOcnTop);
    M.vfm(M.iOcnBot) = (4/3)*pi*(M.r(M.iOcnBot+1)^3-M.rSil^3)/M.V_s(M.iOcnBot);    
    
    % Restore correct temperatures now that energy is back in the right place
    M.T(M.r>=M.rSil & M.r<=M.rOcn) = IN.Tm_ocn;
end


%%%%%%%%%%%%%%%%%%
% Reservoir Melting/Freezing
%%%%%%%%%%%%%%%%%%

if M.vRes>0   
    % Calculate heat extracted from reservoir
    % check whether the interface has warmed above Tm
    dEice_top = (M.T(M.iResTop+1)-M.Tm_res)*M.rho(M.iResTop+1)*M.Cp(M.iResTop+1)*(4/3)*pi*(M.r_s(M.iResTop+1)^3 - M.r_s(M.iResTop)^3);
    dEice_top(dEice_top<0) = 0;  % Ice is allowed to be colder than melting
    if M.T(M.iResTop+1)>M.Tm_res % Ice is not allowed to be warmer
        M.T(M.iResTop+1) = M.Tm_res;
    end

    dEice_bot = (M.T(M.iResBot)-M.Tm_res)*M.rho(M.iResBot)*M.Cp(M.iResBot)*(4/3)*pi*(M.r_s(M.iResBot)^3 - M.r_s(M.iResBot-1)^3);
    dEice_bot(dEice_bot<0) = 0; % Ice is allowed to be colder than melting
    if M.T(M.iResBot)>M.Tm_res  % Ice is not allowed to be warmer
        M.T(M.iResBot) = M.Tm_res;
    end
    
    % Calculate change in melt
    dETop = dE_resTop + dEice_top;
    dEBot = dE_resBot + dEice_bot;
    
    % NOTE: These dvs are for a global spherical shell of melt, NOT for the
    % reservoir! 
    dmTop = dETop/M.LH2O;    % Total change in melt mass
    dvTop = dmTop/M.rhoRes; % Total change in melt volume 
    
    dmBot = dEBot/M.LH2O;    % Total change in melt mass
    dvBot = dmBot/M.rhoRes; % Total change in melt volume
    
    % Change in interface location
    M.rResTop = ((3/(4*pi))*dvTop+M.rResTop^3)^(1/3);
    M.rResBot = (-(3/(4*pi))*dvBot+M.rResBot^3)^(1/3);
    
    % Change in reservoir size and frozen fraction
    M.rRes  = (M.rResTop-M.rResBot)/2; % Reservoir radius
    M.vRes  = (4/3)*pi*M.rRes^3;       % Reservoir volume
    M.fV    = 1-M.rRes^3/IN.rRes^3;    % Reservoir frozen fraction
    M.mRes   = M.vRes*M.rhoRes;        % Reservoir mass
    
    M.iResTop = find((M.rResTop - M.r>0)>0,1,'last'); % Reservoir top interface element index
    M.iResBot = find((M.rResBot - M.r>0)>0,1,'last'); % Reservoir bottom interface element index
    
    if M.rResTop <= M.rResBot 
        % Fully frozen reservoir
        % Energy Accounting
        rho_s = n2sVolumetric(M,M.rho);
        Cp_s  = n2sVolumetric(M,M.Cp);
        
        dE = (4/3)*pi*(M.rResBot^3 - M.rResTop^3)*M.rhoRes*M.LH2O;
        dT = dE/(M.V_s(M.iResTop)*rho_s(M.iResTop)*Cp_s(M.iResTop));
        M.T(M.iResTop) = M.T(M.iResTop)+dT;
        
        % Close out frozen reservoir
        M.vfm(M.r_s>M.rOcn) = 0;
        M.iResBot = M.iResTop;
        rResMean  = mean([M.rResTop,M.rResBot]);
        M.rResTop = rResMean;
        M.rResBot = rResMean;
        M.rRes  = 0;  % Reservoir radius
        M.vRes  = 0;  % Reservoir volume
        M.fV    = 1;  % Reservoir frozen fraction
        M.mRes  = 0; % Reservoir mass
        
    elseif M.rResBot <= M.rOcn
        % Ice shell melted to reservoir depth
        error('Ice Shell Melted Through.');
    else 
        % Set melt fractions
        M.vfm((M.r_s>=M.rResBot) & (M.r_s<=M.rResTop)) = 1;
        M.vfm(M.iResTop) = (4/3)*pi*(M.rResTop^3 - M.r(M.iResTop)^3)/M.V_s(M.iResTop);
        M.vfm(M.iResBot) = 1-(4/3)*pi*(M.rResBot^3 - M.r(M.iResBot)^3)/M.V_s(M.iResBot);
        
        % Restore correct temperatures now that energy is back in the right place
        M.T((M.r>=M.rResBot) & (M.r<=M.rResTop)) = M.Tm_res;
    end    
    
else
    
end



end



















