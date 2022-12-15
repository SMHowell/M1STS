%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IN, OUT, MISC] = outputManager(M,IN,BOD,OUT,MISC)

if IN.outInd * IN.tOut <= M.t
    IN.outInd = IN.outInd + 1;
    
    %%%%%%%%%%%%%%%%%%
    % Save Output
    %%%%%%%%%%%%%%%%%%
    % Note that outputs are initialized in initializeOutputs.m
    OUT.t(IN.outInd)     = M.t;
    OUT.dm_dt(IN.outInd) = M.dm_dt ;
    OUT.Tsurf(IN.outInd) = M.Tsurf;  % Surface Temperature
    OUT.Dice(IN.outInd)  = BOD.R-M.rOcn; % Ice thickness
    
    %%%%%%%%%%%%%%%%%%
    % Display Time
    %%%%%%%%%%%%%%%%%%
    if M.t<IN.kyr2s
        outTime  = M.t/IN.yr2s;
        timeUnit = ' yr';
    elseif M.t<IN.Myr2s
        outTime  = M.t/IN.kyr2s;
        timeUnit = 'kyr';
    elseif M.t<IN.Gyr2s
        outTime  = M.t/IN.Myr2s;
        timeUnit = 'Myr';
    else
        outTime  = M.t/IN.Gyr2s;
        timeUnit = 'Gyr';
    end
    
    
    %%%%%%%%%%%%%%%%%%
    % Time Since Last Output
    %%%%%%%%%%%%%%%%%%
    dt_simu = M.t - M.simuTime_old;
    dt_real = toc - M.realTime_old;
    
    dt_dt   = dt_simu/dt_real;
    if dt_dt<IN.kyr2s
        outRate  = dt_dt/IN.yr2s;
        rateUnit = ' yr/s';
    elseif dt_dt<IN.Myr2s
        outRate  = dt_dt/IN.kyr2s;
        rateUnit = 'kyr/s';
    elseif dt_dt<IN.Gyr2s
        outRate  = dt_dt/IN.Myr2s;
        rateUnit = 'Myr/s';
    else
        outRate  = dt_dt/IN.Gyr2s;
        rateUnit = 'Gyr/s';
    end
    
    M.simuTime_old = M.t;
    M.realTime_old = toc;
    
    %%%%%%%%%%%%%%%%%%
    % Build Output Message
    %%%%%%%%%%%%%%%%%%
    str{1} = [num2str(outTime,'%05.1f') ' ' timeUnit ' Elapsed'];
    str{end+1} = [' @ ' num2str(outRate,'%05.1f') ' ' rateUnit '.'];
    str{end+1} = [' Ice Shell ' num2str(OUT.Dice(IN.outInd)/1e3,'%03.1f') ' km.' ];
    str{end+1} = [' diffMax ' num2str(max(M.fmDiff_H2O),'%01.2f') '.'];
    % str{end+1} = [' T_surf = ' num2str(M.Tsurf,'%05.1f.')];
    str{end+1} = [' ' num2str(100*IN.outInd/IN.NOut,'%04.1f') '% Complete.'];
    disp([str{:}]);
    
    
    %%%%%%%%%%%%%%%%%%
    % Plots
    %%%%%%%%%%%%%%%%%%
    MISC.outTime = outTime;
    MISC.unit = timeUnit;
    
    [MISC] = polarPlotInt(IN,M,OUT,MISC);
    
    if (IN.outInd == IN.NOut)
        save(['Output/' IN.outputName '.mat'],'OUT');
    end
    
end

end























