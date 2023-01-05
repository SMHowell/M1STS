%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = meltingFreezing(M,IN,COMP)


%%%%%%%%%%%%%%%%%%
% Ocean Melting/Freezing
%%%%%%%%%%%%%%%%%%

% Calculate heat extracted from ocean
dE_ocn = (M.T(M.iOcnTop)-M.Tm_ocn)*M.rho(M.iOcnTop)*M.Cp(M.iOcnTop)*(4/3)*pi*(M.r_s(M.iOcnTop)^3 - M.r_s(M.iOcnTop-1)^3);

% Also, check whether the interface has warmed above Tm
dE_ice = (M.T(M.iOcnTop+1)-M.Tm_ocn)*M.rho(M.iOcnTop+1)*M.Cp(M.iOcnTop+1)*(4/3)*pi*(M.r_s(M.iOcnTop+1)^3 - M.r_s(M.iOcnTop)^3);
dE_ice(dE_ice<0) = 0;         % Ice is allowed to be colder than melting
if M.T(M.iOcnTop+1)>M.Tm_ocn % Ice is not allowed to be warmer
    M.T(M.iOcnTop+1) = M.Tm_ocn;
end

% Calculate change in melt
dE = dE_ocn + dE_ice;
dm = dE/IN.L;    % Total change in melt mass
dv = dm/M.rhoOcn; % Total change in melt volume

M.dm_dt = dm/M.dt;

% New interface location
M.rOcnTop = ((3/(4*pi))*dv+M.rOcnTop^3)^(1/3);
M.iOcnTop = find((M.rOcnTop - M.r>0)>0,1,'last');

% Set melt fractions
M.vfm(M.r_s<M.rOcnTop) = 1;
M.vfm(M.r_s>M.rOcnTop) = 0;
M.vfm(M.iOcnTop)       = (4/3)*pi*(M.rOcnTop^3 - M.r(M.iOcnTop)^3)/M.V_s(M.iOcnTop);

% Restore correct temperatures now that energy is back in the right place
M.T(M.T >M.Tm_ocn) = M.Tm_ocn;
M.T(M.r<=M.rOcnTop) = M.Tm_ocn;

% Calculate change in melt
dE = dE_ocn + dE_ice;
dm = dE/IN.L;    % Total change in melt mass
dv = dm/M.rhoOcn; % Total change in melt volume

M.dm_dt = dm/M.dt;

% New interface location
M.rOcnTop = ((3/(4*pi))*dv+M.rOcnTop^3)^(1/3);
M.iOcnTop = find((M.rOcnTop - M.r>0)>0,1,'last');

% Set melt fractions
M.vfm(M.r_s<M.rOcnTop) = 1;
M.vfm(M.r_s>M.rOcnTop) = 0;
M.vfm(M.iOcnTop) = (4/3)*pi*(M.rOcnTop^3 - M.r(M.iOcnTop)^3)/M.V_s(M.iOcnTop);

% Restore correct temperatures now that energy is back in the right place
M.T(M.r<=M.rOcnTop) = M.Tm_ocn;


%%%%%%%%%%%%%%%%%%
% Reservoir Melting/Freezing
%%%%%%%%%%%%%%%%%%
% NOTE: This works by getting the energy extracted from an annulus of water
% spanning the globe. We only use this energy value to calculate change in
% reservoir interface position. By calculating change in reservoir radius
% at the top and bottom, we can estimate overall change in reservoir volume
% effectively.

if M.vRes>0  

    % Grab reservoir properties from composition data
    [T, rho, vlf, vif] = coolFreeze(M, COMP, IN);
    M.vlf = vlf;    % volumic liquid fraction
    M.vif = vif;    % volumic ice fraction (after expansion! vif=V_ice/V_init)

%       Everything commented hereafter is not necessary anymore (I think??)

%     % Calculate heat extracted from reservoir
%     if M.iResTop==M.iResBot
%         % Handle case when reservoir is collapsing
%         dE_resTop = 0.5*(M.T(M.iResTop)-M.Tm_res)*M.rho(M.iResTop)*M.Cp(M.iResTop)*(4/3)*pi*(M.r_s(M.iResTop)^3 - M.r_s(M.iResTop-1)^3);
%         dE_resBot = 0.5*((M.T(M.iResBot+1)-M.Tm_res)*M.rho(M.iResBot+1)*M.Cp(M.iResBot+1)*(4/3)*pi*(M.r_s(M.iResBot+1)^3 - M.r_s(M.iResBot)^3));
%         
%         % check whether the interface has warmed above Tm
%         dE_iceTop = 0.5*(M.T(M.iResTop+1)-M.Tm_res)*M.rho(M.iResTop+1)*M.Cp(M.iResTop+1)*(4/3)*pi*(M.r_s(M.iResTop+1)^3 - M.r_s(M.iResTop)^3);
%         dE_iceTop(dE_iceTop<0) = 0;  % Ice is allowed to be colder than melting
%         if M.T(M.iResTop+1)>M.Tm_res % Ice is not allowed to be warmer
%             M.T(M.iResTop+1) = M.Tm_res;
%         end
%         
%         dE_iceBot = 0.5*(M.T(M.iResBot)-M.Tm_res)*M.rho(M.iResBot)*M.Cp(M.iResBot)*(4/3)*pi*(M.r_s(M.iResBot)^3 - M.r_s(M.iResBot-1)^3);
%         dE_iceBot(dE_iceBot<0) = 0; % Ice is allowed to be colder than melting
%         if M.T(M.iResBot)>M.Tm_res  % Ice is not allowed to be warmer
%             M.T(M.iResBot) = M.Tm_res;
%         end
%         
%     else
%         % Handle nominal case
%         dE_resTop = (M.T(M.iResTop)-M.Tm_res)*M.rho(M.iResTop)*M.Cp(M.iResTop)*(4/3)*pi*(M.r_s(M.iResTop)^3 - M.r_s(M.iResTop-1)^3);
%         dE_resBot = ((M.T(M.iResBot+1)-M.Tm_res)*M.rho(M.iResBot+1)*M.Cp(M.iResBot+1)*(4/3)*pi*(M.r_s(M.iResBot+1)^3 - M.r_s(M.iResBot)^3));
%         
%         % check whether the interface has warmed above Tm
%         dE_iceTop = (M.T(M.iResTop+1)-M.Tm_res)*M.rho(M.iResTop+1)*M.Cp(M.iResTop+1)*(4/3)*pi*(M.r_s(M.iResTop+1)^3 - M.r_s(M.iResTop)^3);
%         dE_iceTop(dE_iceTop<0) = 0;  % Ice is allowed to be colder than melting
%         if M.T(M.iResTop+1)>M.Tm_res % Ice is not allowed to be warmer
%             M.T(M.iResTop+1) = M.Tm_res;
%         end
%         
%         dE_iceBot = (M.T(M.iResBot)-M.Tm_res)*M.rho(M.iResBot)*M.Cp(M.iResBot)*(4/3)*pi*(M.r_s(M.iResBot)^3 - M.r_s(M.iResBot-1)^3);
%         dE_iceBot(dE_iceBot<0) = 0; % Ice is allowed to be colder than melting
%         if M.T(M.iResBot)>M.Tm_res  % Ice is not allowed to be warmer
%             M.T(M.iResBot) = M.Tm_res;
%         end
%     end
%     
%     % Calculate change in melt
%     dETop = dE_resTop + dE_iceTop;
%     dEBot = dE_resBot + dE_iceBot;
%     
%     % NOTE: These dvs are for a global spherical shell of melt, NOT for the
%     % reservoir! 
%     dmTop = M.DeltaEtop/IN.L;       % Total change in melt mass
%     dvTop = dmTop/M.rhoRes;         % Total change in melt volume 
%     
%     dmBot = M.DeltaEbot/IN.L;       % Total change in melt mass
%     dvBot = dmBot/M.rhoRes;         % Total change in melt volume
%     
%     % Change in interface location
%     M.rResTop = ((3/(4*pi))*dvTop+M.rResTop^3)^(1/3);
%     M.rResBot = (-(3/(4*pi))*dvBot+M.rResBot^3)^(1/3);
%     
%     % Change in reservoir size and frozen fraction
%     M.rRes  = (M.rResTop-M.rResBot)/2; % Reservoir radius
%     M.vRes  = (4/3)*pi*M.rRes^3;       % Reservoir volume
%     M.fV    = 1-M.rRes^3/IN.rRes^3;    % Reservoir frozen fraction
%     M.mRes   = M.vRes*M.rhoRes;        % Reservoir mass
%     
%     M.iResTop = find((M.rResTop - M.r>0)>0,1,'last'); % Reservoir top interface element index
%     M.iResBot = find((M.rResBot - M.r>0)>0,1,'last'); % Reservoir bottom interface element index


    % Update parameters
    M.rhoRes    = rho;
    M.fV        = 1 - vlf;
    M.Tm_res    = T;

    % Change in reservoir size
    vResOld     = M.vRes;
    M.vRes      = (4/3)*pi*IN.rRes^3 * vlf;     % New reservoir volume
    vResDiff    = M.vRes - vResOld;
    vDiffTop    = vResDiff * M.ratioEtop;       % Volume diff at res top
    vDiffBot    = vResDiff * M.ratioEbot;       % Volume diff at res bottom (not sure these are very clean...)
    M.mRes      = M.vRes * M.rhoRes;            % New reservoir mass
    M.rRes      = (3/4/pi*M.vRes)^(1/3);        % Reservoir radius
    M.rResTop   = M.rResTop + sign(vDiffTop)*(abs(vDiffTop)*3/4/pi)^(1/3);  % New top interface                    
    M.rResBot   = M.rResBot + sign(vDiffBot)*(abs(vDiffBot)*3/4/pi)^(1/3);  % New bottom interface
    M.iResTop   = find((M.rResTop - M.r>0)>0,1,'last');   % Reservoir top interface element index
    M.iResBot   = find((M.rResBot - M.r>0)>0,1,'last');   % Reservoir bottom interface element index

    
    if M.rResTop <= M.rResBot 
        % Fully frozen reservoir
        % Energy Accounting
        rho_s = n2sVolumetric(M,M.rho);
        Cp_s  = n2sVolumetric(M,M.Cp);
        
        dE = (4/3)*pi*(M.rResBot^3 - M.rResTop^3)*M.rhoRes*IN.L;
        dT = dE/(M.V_s(M.iResTop)*rho_s(M.iResTop)*Cp_s(M.iResTop));
        M.T(M.iResTop) = M.T(M.iResTop)+dT;
        
        % Close out frozen reservoir
        M.vfm(M.r_s>M.rOcnTop) = 0;
        M.iResBot = M.iResTop;
        rResMean  = mean([M.rResTop,M.rResBot]);
        M.rResTop = rResMean;
        M.rResBot = rResMean;
        M.rRes  = 0;  % Reservoir radius
        M.vRes  = 0;  % Reservoir volume
        M.fV    = 1;  % Reservoir frozen fraction
        M.mRes  = 0;  % Reservoir mass
        
    elseif M.rResBot <= M.rOcnTop
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

end
end



















