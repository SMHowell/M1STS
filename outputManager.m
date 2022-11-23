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

if IN.outInd * IN.tOut < M.t
    IN.outInd = IN.outInd + 1;
    
    %%%%%%%%%%%%%%%%%%
    % Save Output
    %%%%%%%%%%%%%%%%%%
    % Note that outputs are initialized in initializeOutputs.m
    OUT.t(IN.outInd)     = M.t;
    OUT.dm_dt(IN.outInd) = M.dm_dt ;
    OUT.Tsurf(IN.outInd) = M.Tsurf;  % Surface Temperature
    OUT.Dice(IN.outInd)  = BOD.R-M.rOcn; % Ice thickness
    
    if M.fV>0 % Only output of theres a reservoir, else zeros
        OUT.vRes(IN.outInd)  = M.vRes;   % Reservoir volume
        OUT.rRes(IN.outInd)  = M.rRes;   % Reservoir radius
        OUT.zRes(IN.outInd)  = BOD.R-M.rResTop; % Reservoir Depth
        OUT.fV(IN.outInd)    = M.fV;     % Reservoir frozen fraction
        OUT.Tmelt(IN.outInd) = M.Tm_res; % Reservoir melting temperature
    end
    
    %%%%%%%%%%%%%%%%%%
    % Display Time
    %%%%%%%%%%%%%%%%%%
    if M.t<IN.kyr2s
        outTime = M.t/IN.yr2s;
        unit    = ' yr';
    elseif M.t<IN.Myr2s
        outTime = M.t/IN.kyr2s;
        unit    = 'kyr';
    elseif M.t<IN.Gyr2s
        outTime = M.t/IN.Myr2s;
        unit    = 'Myr';
    else
        outTime = M.t/IN.Gyr2s;
        unit    = 'Gyr';
    end
    
    %%%%%%%%%%%%%%%%%%
    % Build Output Message
    %%%%%%%%%%%%%%%%%%
    str{1} = [num2str(outTime,'%05.1f') ' ' unit ' Elapsed.'];
    str{end+1} = [' T_surf = ' num2str(M.Tsurf,'%05.1f.')];
    str{end+1} = [' ' num2str(100*IN.outInd/IN.NOut,'%04.1f') '% Complete.'];
    disp([str{:}]);
    
    %%%%%%%%%%%%%%%%%%
    % Plots
    %%%%%%%%%%%%%%%%%%
    MISC.outTime = outTime;
    MISC.unit = unit;
    [MISC] = polarPlot(IN,M,OUT,MISC);
    
    
    if (IN.outInd == IN.NOut)
        if IN.resOn
            outputName = ['H0' num2str(IN.H0_ice/1e3,'%0.1f') '_rRes' num2str(IN.rRes/1e3,'%0.1f') '_zRes' num2str(IN.zResTop/1e3,'%0.1f')];
            save(['Output/' outputName '.mat'],'OUT');
        else
            outputName = ['H0' num2str(IN.H0_ice/1e3,'%0.1f') '_resOff'];
            save(['Output/' outputName '.mat'],'OUT');
        end
    end
    
end

end























