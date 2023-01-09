%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function COMP = getCompData(IN)

% reads the name of sheets in composition data file
COMP.sheets = sheetnames('CompData.xlsx');

i = IN.simu;

% possible to make iterative to use all compositions:
% for i = 1:length(COMP.sheets)

    % reads composition data sheets
    COMP.Table{i}   = readtable('CompData.xlsx','sheet', COMP.sheets{i}, 'PreserveVariableNames',true);

    % starts at liquidus:
    start = 1;
    
    % temperature cell array
    COMP.temps{i}   = table2array(COMP.Table{i}(start:end,1)).';
    
    % liquid volume cell array
    COMP.V_l{i}     = table2array(COMP.Table{i}(start:end,2)).'/1e6;
    
    % liquid density
    COMP.rho_l{i}   = table2array(COMP.Table{i}(start:end,3)).'*1000;
    
    % ice volume
    COMP.V_ice{i}   = table2array(COMP.Table{i}(start:end,4)).'/1e6;
    
    % making volumes cumulative for the Fractional Freezing cases
    if contains(COMP.sheets{i}, 'Frac') == 1
        COMP.V_ice{i} = cumsum(COMP.V_ice{i});
    end
    
    % volumic liquid fraction
    % COMP.V_init{i}  = table2array(COMP.Table{i}(1,2)) / 1e6;
    COMP.V_init{i}  = 1000 / 1e6;   % Ask Marc and Mariam if the initial volume is always 1000 cm^3 ???
    COMP.vlf{i}     = COMP.V_l{i} ./ COMP.V_init{i}; % issue with datasets where V_i is well bellow 1e-3!!!
    
    if isempty(COMP.vlf{i}(COMP.vlf{i}>1)) == false
        COMP.vlf{i}(COMP.vlf{i}>1) = 1; % "noise" in comp. data
    end
    
    % volumic frozen fraction before expansion
    COMP.vff{i}     = 1 - COMP.vlf{i}; 

    % volumic ice fraction after expansion
    COMP.vif{i}     = COMP.V_ice{i} ./ COMP.V_init{i}; % issue with datasets where V_i is well bellow 1e-3!!!
    
    % number of rows
    nrows           = length(COMP.temps{i});
    
    % reads salt concentrations
    COMP            = getSalts(COMP, i, nrows);

% end

disp('Composition data imported')

end
