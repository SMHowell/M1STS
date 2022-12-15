%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = meltingFreezing(M,IN,BOD)

%%%%%%%%%%%%%%%%%%
% Ocean Melting/Freezing
%%%%%%%%%%%%%%%%%%
% Only execute if there is a hydrosphere
if M.vH2O > 0
    
    %%%%%%%
    % Ocean radiogenic heating
    %%%%%%%
    % Radiate heat in the ocean from leeched 40K
    dE_rad = M.dt * M.Prad_ocn;
    
    
    %%%%%%%
    % Handle energy sources/sinks to ocean.
    %%%%%%%
    if M.iOcnTop==M.iOcnBot
        % Consider the case where the ocean is less than one element deep or
        % fully frozen
        if (M.vOcn==0)
            % No ocean. See if one should form. In this first formation step,
            % we are only concerned with whether the seafloor might be warm
            % enough to begin melting ice.
            
            % Apply ocean radiogenic heat to seafloor
            ind = M.iOcnBot;
            dT  = dE_rad/(2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
                (M.dr(ind-1)+M.dr(ind)));
            M.T(M.iOcnBot) = M.T(M.iOcnBot) + dT;
            
            % Radio energy is now used up
            dE_rad = 0;
            
            % Check if seafloor should cause melting in direct contact with ice.
            ind = M.iOcnBot;
            dE  = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
                (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-IN.Tm_ocn);
            dE(dE<0) = 0; % Only care if this is positive
            
            % Because we withdrew the extra energy, make sure to set
            % temperature accordingly.
            if dE>0
                M.T(ind) = IN.Tm_ocn;
            end
            
        else
            % Ocean exists, but is less than 1 element thick
            % All the seafloor heat beyond that required to reach the melting
            % temperature goes into the ocean, and is used to melt.
            ind = M.iOcnBot;
            dE  = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
                (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-IN.Tm_ocn);
            
            % Because we withdrew the extra energy, make sure to set
            % temperature accordingly.
            M.T(ind) = IN.Tm_ocn;
            
        end
    else
        % Consider the case where the ocean is more than one element thick.
        
        % All the seafloor heat beyond that required to reach the melting
        % temperature goes into the ocean, and is used to melt.
        ind   = M.iOcnBot;
        dE_in = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
            (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-IN.Tm_ocn);
        
        % Now, account for the heat lost from the ocean to the ice shell, which
        % will appear as cooling in the fully-ocean element beneath the
        % interface.
        ind    = M.iOcnTop;
        dE_out = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
            (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-IN.Tm_ocn);
        
        % Because we withdrew the extra energy, make sure to set
        % temperature accordingly.
        M.T(M.iOcnBot:M.iOcnTop) = IN.Tm_ocn;
        
        % Sum energy change
        dE = dE_in + dE_out;
    end
    % Add any remaining radiative energy
    dE = dE + dE_rad;
    
    
    %%%%%%%
    % Interstitial melting
    %%%%%%%
    % Now check for any melting within the icy shell (e.g. from tidal heating).
    % If so, extract it to the ocean.
    
    % Accumulate excess energy from ice warmer than melting temp
    ind = find(M.r>=M.rOcn & M.r<=BOD.R);
    dE_in = 2*pi*M.rho(ind).*M.Cp(ind).*M.r(ind).^2 ...
        .*(M.dr(ind-1)+M.dr(ind)).*(M.T(ind)-IN.Tm_ocn);
    dE_in(dE_in<0) = 0;
    
    % Add to running total
    dE = dE + sum(dE_in);
    
    % Because we withdrew the extra energy, make sure to set
    % temperature accordingly.
    T = M.T(ind);
    T(dE_in>0) = IN.Tm_ocn;
    M.T(ind)   = T;
    
    
    %%%%%%%
    % Melt / Freeze
    %%%%%%%
    if abs(dE)>0 % On
        % Get interface location and energy density of phase change, U [J/m^3]
        [U,rOcn_new] = getOceanRadius(M,BOD,IN,dE);
        
        % Get new ocean volume
        vOcn_new = (4/3)*pi*(rOcn_new^3 - M.rSil^3);
        
        % Total change in melt mass
        dm = (vOcn_new - M.vOcn)*MAT.H2O.m.rho0;
        
        % Melting Rate
        M.dm_dt = dm/M.dt;
        
        % Update ocean geometry
        M.rOcn = rOcn_new;
        M.vOcn = vOcn_new;
        M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last');
        
        % Set melt fractions on elements
        M.vfm(M.r_s>M.rSil & M.r_s<M.rOcn) = 1;
        M.vfm(M.r_s>M.rOcn | M.r_s<M.rSil) = 0;
        M.vfm(M.iOcnTop) = (4/3)*pi*(M.rOcn^3 - M.r(M.iOcnTop)^3)/M.V_s(M.iOcnTop);
        M.vfm(M.iOcnBot) = min((4/3)*pi*(M.rOcn^3 - M.rSil^3)/M.V_s(M.iOcnBot),1);
        
        % Check to see if ocean has frozen out completely
        if M.rOcn<M.rSil
            % Find excess energy that went into melting
            dE = (4/3)*pi*(M.rOcn^3 - M.rSil^3)*U;
            
            % Instead use to cool seafloor
            ind = M.iOcnBot;
            dT  = dE/(2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
                (M.dr(ind-1)+M.dr(ind)));
            M.T(ind) = M.T(ind) + dT;
            
            % Now that energy was spent appropraitely, fix depth
            M.rOcn = M.rSil;
            M.vOcn = 0;
            M.iOcnTop = M.iOcnBot;
        end
    else
        M.dm_dt = 0;
    end
    
end % end if hydrosphere
end



















