%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NOT FUNCTIONNAL YET


function M = Eruption(COMP,M)

if M.vRes>0

    % freezing-induced overpressure in the reservoir:
    V_ice           = M.vRes_init * M.vif;       % volume of ice (after expansion)
    Vl_compressed   = M.vRes_init - V_ice;       % volume of liquid once compressed (after ice expansion)
    Vl              = M.vRes_init * (1-M.vlf);   % volume of liquid if not pressurized
    M.deltaP        = -1/IN.X * log(Vl_compressed/Vl);
    
    % check if it overcomes threshold overpressure:
    if deltaP >= ICE.deltaP_c
        M.eruption = 1;
        OUT.V_erupt = OUT.Vl - OUT.Vl_compressed;
    
        fprintf('Erupted volume: %f m^3 \n',OUT.V_erupt);
    
        % find index of closest remaining energy fraction from composition data
        i_low = find(COMP.fE_rmn{IN.simu} > OUT.fE_rmn, 1, 'last');
        
        % interpolate values of temperature, volumic frozen fraction and density
        % associated
        if i_low == 1 
            indexInterp = [i_low, i_low+1, i_low+2];
        else
            indexInterp = [i_low-1, i_low, i_low+1, i_low+2];
        end
        % Erupted composition:
        OUT.Ca  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Ca{IN.simu}(indexInterp),OUT.fE_rmn,'spline');
        OUT.Mg  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Mg{IN.simu}(indexInterp),OUT.fE_rmn,'spline');
        OUT.Na  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Na{IN.simu}(indexInterp),OUT.fE_rmn,'spline');
        OUT.K   = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.K{IN.simu}(indexInterp),OUT.fE_rmn,'spline');
        OUT.Cl  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Cl{IN.simu}(indexInterp),OUT.fE_rmn,'spline');
        OUT.S   = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.S{IN.simu}(indexInterp),OUT.fE_rmn,'spline');
        OUT.C   = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.C{IN.simu}(indexInterp),OUT.fE_rmn,'spline');
        OUT.Si  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Si{IN.simu}(indexInterp),OUT.fE_rmn,'spline');
    
        fprintf('Erupted composition [mol/kg of water]: %f Ca + %f Mg + %f Na + %f K + %f Cl + %f S + %f C + %f Si\n',OUT.Ca, OUT.Mg, OUT.Na, OUT.K, OUT.Cl, OUT.S, OUT.C, OUT.Si);
    
    end
end

end