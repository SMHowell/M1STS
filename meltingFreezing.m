%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = meltingFreezing(M,IN,BOD,MAT)

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
                (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-MAT.H2O.Tm0);
            dE(dE<0) = 0; % Only care if this is positive
            
            % Because we withdrew the extra energy, make sure to set
            % temperature accordingly.
            if dE>0
                M.T(ind) = MAT.H2O.Tm0;
            end
            
        else
            % Ocean exists, but is less than 1 element thick
            % All the seafloor heat beyond that required to reach the melting
            % temperature goes into the ocean, and is used to melt.
            ind = M.iOcnBot;
            dE  = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
                (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-MAT.H2O.Tm0);
            
            % Because we withdrew the extra energy, make sure to set
            % temperature accordingly.
            M.T(ind) = MAT.H2O.Tm0;
            
        end
    else
        % Consider the case where the ocean is more than one element thick.
        
        % All the seafloor heat beyond that required to reach the melting
        % temperature goes into the ocean, and is used to melt.
        ind   = M.iOcnBot;
        dE_in = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
            (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-MAT.H2O.Tm0);
        
        % Now, account for the heat lost from the ocean to the ice shell, which
        % will appear as cooling in the fully-ocean element beneath the
        % interface.
        ind    = M.iOcnTop;
        dE_out = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
            (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-MAT.H2O.Tm0);
        
        % Because we withdrew the extra energy, make sure to set
        % temperature accordingly.
        M.T(M.iOcnBot:M.iOcnTop) = MAT.H2O.Tm0;
        M.T(end-1:end) = M.Tsurf;
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
        .*(M.dr(ind-1)+M.dr(ind)).*(M.T(ind)-MAT.H2O.Tm0);
    dE_in(dE_in<0) = 0;
    
    % Add to running total
    dE = dE + sum(dE_in);
    
    % Because we withdrew the extra energy, make sure to set
    % temperature accordingly.
    T = M.T(ind);
    T(dE_in>0) = MAT.H2O.Tm0;
    M.T(ind)   = T;
    
    
    %%%%%%%
    % Melt / Freeze
    %%%%%%%
    if abs(dE)>0 % On
        % Get interface location and energy density of phase change, U [J/m^3]
        [U,rOcn_new] = getOceanRadius(M,BOD,IN,dE,MAT);
        
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
        vfm = M.mat.fm_s(M.mat.iH2Omelt,:);
        vfm(M.r_s>M.rSil & M.r_s<M.rOcn) = 1;
        vfm(M.r_s>M.rOcn | M.r_s<M.rSil) = 0;
        vfm(M.iOcnTop) = (4/3)*pi*(M.rOcn^3 - M.r(M.iOcnTop)^3)/M.V_s(M.iOcnTop);
        vfm(M.iOcnBot) = min((4/3)*pi*(M.rOcn^3 - M.rSil^3)/M.V_s(M.iOcnBot),1);
        dvfm =  vfm-M.mat.fm_s(M.mat.iH2Omelt,:);
        
        M.mat.fm_s(M.mat.iH2Omelt,:)  = M.mat.fm_s(M.mat.iH2Omelt,:)  + dvfm;
        M.mat.fm_s(M.mat.iH2Osolid,:) = M.mat.fm_s(M.mat.iH2Osolid,:) - dvfm;
        
        rhoFull_s = n2sVolumetric(M,M.mat.rhoFull);
        drho_s  = sum(M.mat.fm_s./rhoFull_s).^-1-M.rho_s;
        drho    = s2nVolumetric(M,drho_s);
        M.rho   = M.rho+drho;
        M.rho_s = M.rho_s+drho_s;
        M.mat.fV_s = (M.rho_s .* M.mat.fm_s ./ rhoFull_s); 
        
        % Make sure nothing got weird 
        M.mat.fm_s(M.mat.fm_s<0) = 0;
        M.mat.fm_s(M.mat.fV_s<0) = 0;
        M.mat.fm_s(M.mat.fm_s>1) = 1;
        M.mat.fm_s(M.mat.fV_s>1) = 1;
        
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
        if M.rOcn>=BOD.R-IN.Hmin_H2O 
            % Dont let thin below critical value. Return heatflow to
            % surface temperature estimation
            dV = (4/3)*pi*(M.rOcn^3 - (BOD.R-IN.Hmin_H2O)^3);
            dm = dV .* MAT.H2O.m.rho0;
            dE = dm * MAT.H2O.L;
            
            % Now that energy was spent appropraitely, fix depth and temp
            M.rOcn = BOD.R-IN.Hmin_H2O;
            M.T(M.iOcnBot+1:M.iOcnTop) = MAT.H2O.Tm0;
            M.Tsurf = MAT.H2O.Tm0;
            
            % Fix surface temp. based on imposing minimum thickness.
            qSurfAddtl = (dE/M.dt)/BOD.Asurf;
            M = getMeltFixHeatFlux(M,BOD,IN,MAT,qSurfAddtl);
            
        end
    else
        M.dm_dt = 0;
    end
    
end % end if hydrosphere





% %%%%%%%%%%%%%%%%%%
% % Iron Melting/Freezing
% %%%%%%%%%%%%%%%%%%
% % Only execute if there is a core
% if M.vIrn > 0
%     
%     
%     %%%%%%%
%     % Handle energy sources/sinks to outer core
%     %%%%%%%
%     if M.iIocTop==M.iIocBot
%         % Consider the case where the OC is less than one element deep or
%         % fully frozen
%         if (M.vIoc==0)
%             % No core. See if one should form. In this first formation step,
%             % we are only concerned with whether the core might be warm
%             % enough to begin melting iron.
%             ind = M.iIocBot;
%             if ind==1
%                 dE  = 2*pi*M.rho(ind)*M.Cp(ind)*M.rIrn^2*...
%                     (M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%             else
%                 dE  = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
%                     (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%             end
%             dE(dE<0) = 0; % Only care if this is positive
%             
%             % Because we withdrew the extra energy, make sure to set
%             % temperature accordingly.
%             if dE>0
%                 M.T(ind) = MAT.IRN.Tm0;
%             end
%             
%         else
%             % Ocean exists, but is less than 1 element thick
%             % All the seafloor heat beyond that required to reach the melting
%             % temperature goes into the ocean, and is used to melt.
%             ind = M.iIocBot;
%             if ind==1
%                 dE  = 2*pi*M.rho(ind)*M.Cp(ind)*M.rIrn^2*...
%                     (M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%             else
%                 dE  = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
%                     (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%             end
%             % Because we withdrew the extra energy, make sure to set
%             % temperature accordingly.
%             M.T(ind) = MAT.IRN.Tm0;
%         end
%     else
%         % Consider the case where the ocean is more than one element thick.
%         
%         % All the seafloor heat beyond that required to reach the melting
%         % temperature goes into the ocean, and is used to melt.
%         ind   = M.iIocBot;
%         if ind==1
%             dE_in  = 2*pi*M.rho(ind)*M.Cp(ind)*M.rIrn^2*...
%                 (M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%         else
%             dE_in  = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
%                 (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%         end
%         
%         % Now, account for the heat lost from the ocean to the ice shell, which
%         % will appear as cooling in the fully-ocean element beneath the
%         % interface.
%         ind    = M.iIocTop;
%         if ind==1
%             dE_out  = 2*pi*M.rho(ind)*M.Cp(ind)*M.rIrn^2*...
%                 (M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%         else
%             dE_out  = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
%                 (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%         end
%         
%         % Because we withdrew the extra energy, make sure to set
%         % temperature accordingly.
%         M.T(M.iIocBot:M.iIocTop) = MAT.IRN.Tm0;
%         
%         % Sum energy change
%         dE = dE_in + dE_out;
%     end  
%     
%     %%%%%%%
%     % Interstitial melting
%     %%%%%%%
%     % Now check for any melting within the core (e.g. from tidal heating).
%     % If so, extract it to the OC.
%     
%     % Accumulate excess energy from ice warmer than melting temp
%     ind = find(M.r<=M.rIoc);
%     if ind==1
%         dE_in  = 2*pi*M.rho(ind)*M.Cp(ind)*M.rIrn^2*...
%             (M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%     else
%         dE_in  = 2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
%             (M.dr(ind-1)+M.dr(ind))*(M.T(ind)-MAT.IRN.Tm0);
%     end
%     dE_in(dE_in<0) = 0;
%     
%     % Add to running total
%     dE = dE + sum(dE_in);
%     
%     % Because we withdrew the extra energy, make sure to set
%     % temperature accordingly.
%     T = M.T(ind);
%     T(dE_in>0) = MAT.IRN.Tm0;
%     M.T(ind)   = T;
%     
%     
%     %%%%%%%
%     % Melt / Freeze
%     %%%%%%%
%     if abs(dE)>0 % On
%         % Get interface location and energy density of phase change, U [J/m^3]
%         [U,rIoc_new] = getCoreRadius(M,BOD,IN,dE,MAT);
%         
%         % Get new ocean volume
%         vIoc_new = (4/3)*pi*(M.rIrn^3 - rIoc_new^3);
%         
%         % Total change in melt mass
%         dm = (vIoc_new - M.vIoc)*MAT.IRN.m.rho0;
%         
%         % Melting Rate
%         M.dm_dt = dm/M.dt;
%         
%         % Update ocean geometry
%         M.rIoc = rIoc_new;
%         M.vIoc = vIoc_new;
%         M.iIocTop = find((M.rIoc - M.r)>=0,1,'last');
%         
%         % Set melt fractions on elements
%         vfm = M.mat.fm_s(M.mat.iIrnMelt,:);
%         vfm(M.r_s>M.rIoc & M.r_s<M.rIrn) = 1;
%         vfm(M.r_s>M.rIrn | M.r_s<M.rIoc) = 0;
%         vfm(M.iIocTop) = (4/3)*pi*(M.rIrn^3-M.r(M.iIocTop)^3)/M.V_s(M.iIocTop);
%         vfm(M.iIocBot) = min((4/3)*pi*(M.r(M.iIocBot+1)^3-M.rIoc^3)/M.V_s(M.iIocBot),1);
%         dvfm =  vfm - M.mat.fm_s(M.mat.iIrnMelt,:);
%         
%         M.mat.fm_s(M.mat.iIrnMelt,:)  = M.mat.fm_s(M.mat.iIrnMelt,:)  + dvfm;
%         M.mat.fm_s(M.mat.iH2Osolid,:) = M.mat.fm_s(M.mat.iH2Osolid,:) - dvfm;
%         
%         rhoFull_s = n2sVolumetric(M,M.mat.rhoFull);
%         drho_s  = sum(M.mat.fm_s./rhoFull_s).^-1-M.rho_s;
%         drho    = s2nVolumetric(M,drho_s);
%         M.rho   = M.rho+drho;
%         M.rho_s = M.rho_s+drho_s;
%         M.mat.fV_s = (M.rho_s .* M.mat.fm_s ./ rhoFull_s); 
%         
%         % Make sure nothing got weird 
%         M.mat.fm_s(M.mat.fm_s<0) = 0;
%         M.mat.fm_s(M.mat.fV_s<0) = 0;
%         M.mat.fm_s(M.mat.fm_s>1) = 1;
%         M.mat.fm_s(M.mat.fV_s>1) = 1;
%         
%         % Check to see if outer core has frozen out completely
%         if M.rIoc<M.rIrn
%             % Find excess energy that went into melting
%             dE = (4/3)*pi*(M.rIoc^3 - M.rIrn^3)*U;
%             
%             % Instead use to cool seafloor
%             ind = M.iIocBot;
%             if ind==1
%                 dT  = dE/(2*pi*M.rho(ind)*M.Cp(ind)*M.rIrn^2*...
%                 (M.dr(ind)));
%             else
%                 dT  = dE/(2*pi*M.rho(ind)*M.Cp(ind)*M.r(ind)^2*...
%                     (M.dr(ind-1)+M.dr(ind)));
%             end
%             M.T(ind) = M.T(ind) + dT;
%             
%             % Now that energy was spent appropraitely, fix depth
%             M.rIoc = M.rIrn;
%             M.vIoc = 0;
%             M.iIocTop = M.iIocBot;
%         end
%     end
%     
% end % end if core

end















