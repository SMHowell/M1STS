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


function M = Eruption(IN,COMP,M)

if M.vRes>0

    % freezing-induced overpressure in the reservoir:

    % Vl = M.vRes_init * M.vlf;
    % vlfRef = Vl / M.vRef;
    % vffRef = 1 - vlfRef;

    M.ViceTot         = M.vRes_init * M.vif;       % volume of ice (after expansion)
    M.Vice            = M.ViceTot - M.Vice_old;    % volume of ice since last eruption
    M.Vl_compressed   = M.vRes_old - M.Vice;         % volume of liquid once compressed (after ice expansion)
    M.Vl              = M.vRes_old - M.vRes_init*(1-M.vlf);       % volume of liquid if not pressurized
    M.deltaP          = -1/IN.X * log(M.Vl_compressed/M.Vl);
    
    % check if it overcomes threshold overpressure:
    if M.deltaP >= M.DeltaPc
        M.eruption  = M.eruption + 1;
        M.V_erupt   = M.Vl - M.Vl_compressed;
    
        fprintf('Erupted volume: %f km^3 \n',M.V_erupt*1e-9);
        fprintf('Freezing time : %f yr \n',(M.t-IN.tRes)/365.25/24/3600);
    
        % find index of closest remaining energy fraction from composition data
        i_low = find(COMP.fE_rmn{IN.simu} >= M.fE_rmn, 1, 'last');
        i_max = length(COMP.fE_rmn{IN.simu});
        
        % interpolate values of temperature, volumic frozen fraction and density
        % associated
        if i_low == 1 
        indexInterp = [i_low, i_low+1];
        elseif i_low+1 == i_max
            indexInterp = [i_low-1, i_low, i_low+1];
        elseif i_low == i_max
            indexInterp = [i_low-2, i_low-1, i_low];
        else
            indexInterp = [i_low-1, i_low, i_low+1, i_low+2];
        end

        % Erupted composition:
        M.Ca  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Ca{IN.simu}(indexInterp),M.fE_rmn,'spline');
        M.Mg  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Mg{IN.simu}(indexInterp),M.fE_rmn,'spline');
        M.Na  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Na{IN.simu}(indexInterp),M.fE_rmn,'spline');
        M.K   = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.K{IN.simu}(indexInterp),M.fE_rmn,'spline');
        M.Cl  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Cl{IN.simu}(indexInterp),M.fE_rmn,'spline');
        M.S   = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.S{IN.simu}(indexInterp),M.fE_rmn,'spline');
        M.C   = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.C{IN.simu}(indexInterp),M.fE_rmn,'spline');
        M.Si  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Si{IN.simu}(indexInterp),M.fE_rmn,'spline');
    
        fprintf('Input composition   [mol/kg of water]: %f Ca + %f Mg + %f Na + %f K + %f Cl + %f S + %f C + %f Si\n',COMP.Ca{IN.simu}(1), COMP.Mg{IN.simu}(1), COMP.Na{IN.simu}(1), COMP.K{IN.simu}(1), COMP.Cl{IN.simu}(1), COMP.S{IN.simu}(1), COMP.C{IN.simu}(1), COMP.Si{IN.simu}(1));
        fprintf('Erupted composition [mol/kg of water]: %f Ca + %f Mg + %f Na + %f K + %f Cl + %f S + %f C + %f Si\n',M.Ca, M.Mg, M.Na, M.K, M.Cl, M.S, M.C, M.Si);
        
        if M.eruption > 50
            error('50 eruptions')
        end

        % update remaining liquid volume:
        M.Vice_old = M.ViceTot;
        M.vRes_old = M.Vl_compressed;


    end
end

end