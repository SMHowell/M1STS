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


function [M, OUT] = Eruption(IN,COMP,M,OUT)

if M.vRes>0 && M.canErupt

    % freezing-induced overpressure in the reservoir    
    M.ViceTot         = M.vRes_init * M.vif;                    % total volume of ice (after expansion)
    M.Vice            = M.ViceTot - M.Vice_old;                 % volume of ice since last eruption
    M.Vl_compressed   = M.vRes_old - M.Vice;                    % volume of liquid once compressed (after ice expansion)
    M.Vl              = M.vRes_old - M.vRes_init*(1-M.vlf);     % volume of liquid if not pressurized
    M.deltaP          = -1/IN.X * log(M.Vl_compressed/M.Vl);    % overpressure
    
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
    M.compLabels = {'Ca','Mg','Na','K','Cl','S','C','Si'};
    M.compCa  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Ca{IN.simu}(indexInterp),M.fE_rmn,'spline');
    M.compMg  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Mg{IN.simu}(indexInterp),M.fE_rmn,'spline');
    M.compNa  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Na{IN.simu}(indexInterp),M.fE_rmn,'spline');
    M.compK   = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.K{IN.simu}(indexInterp),M.fE_rmn,'spline');
    M.compCl  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Cl{IN.simu}(indexInterp),M.fE_rmn,'spline');
    M.compS   = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.S{IN.simu}(indexInterp),M.fE_rmn,'spline');
    M.compC   = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.C{IN.simu}(indexInterp),M.fE_rmn,'spline');
    M.compSi  = interp1(COMP.fE_rmn{IN.simu}(indexInterp),COMP.Si{IN.simu}(indexInterp),M.fE_rmn,'spline');

    % check if it overcomes threshold overpressure:
    if M.deltaP >= M.DeltaPc

        M.V_erupt   = M.Vl - M.Vl_compressed;
    
        fprintf('Erupted volume: %f km^3 \n',M.V_erupt*1e-9);
        fprintf('Freezing time : %f yr \n',(M.t-IN.tRes)/365.25/24/3600);
        fprintf('Input composition   [mol/kg of water]: %f Ca + %f Mg + %f Na + %f K + %f Cl + %f S + %f C + %f Si\n',COMP.Ca{IN.simu}(1), COMP.Mg{IN.simu}(1), COMP.Na{IN.simu}(1), COMP.K{IN.simu}(1), COMP.Cl{IN.simu}(1), COMP.S{IN.simu}(1), COMP.C{IN.simu}(1), COMP.Si{IN.simu}(1));
        fprintf('Erupted composition [mol/kg of water]: %f Ca + %f Mg + %f Na + %f K + %f Cl + %f S + %f C + %f Si\n',M.compCa, M.compMg, M.compNa, M.compK, M.compCl, M.compS, M.compC, M.compSi);

        if M.eruption > IN.nErupt
            error('Too many eruptions')
        end

        % update remaining liquid volume:
        M.Vice_old = M.ViceTot;
        M.vRes_old = M.Vl_compressed;

        OUT.eruptTimes(M.eruption) = M.t;
        OUT.eruptV(M.eruption)     =  M.Vice;
        
        M.eruption  = M.eruption + 1;

    end
end

end