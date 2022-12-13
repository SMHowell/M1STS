%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elodie Lesage & Samuel M. Howell & Julia W. Miller
% NASA Jet Propulsion Laboratory
% (c)2022 California Institute of Technology
% All rights reserved.
% elodie.lesage@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function COMP = getCompData(RES, ICE, COMP)

% reads the name of sheets in composition data file
COMP.sheets = sheetnames('CompData.xlsx');

i = COMP.simu;

% for i = 1:length(COMP.sheets)

    % reads composition data sheets
    COMP.Table{i}   = readtable('CompData.xlsx','sheet', COMP.sheets{i}, 'PreserveVariableNames',true);
    
    % starts at given T
    % start = find(table2array(COMP.Table{i}(:,1)) <= 273, 1, 'first');
    % starts at liquidus
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
    COMP.V_init{i}  = table2array(COMP.Table{i}(1,2)) / 1e6;
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the energy loss from liquid at 
    % each temperature step and how it's spent
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % total energy in liquid at liquidus
    E_init              = COMP.rho_l{i}(1)*COMP.V_init{i}(1) * (COMP.temps{i}(1)*RES.cp + ICE.L);
    
    prevT               = [];
    prevT(1)            = 0;                            % temperature at previous temp. step
    prevT(2:nrows)      = COMP.temps{i}(1:nrows-1);          
    prevV_l             = [];
    prevV_l(1)          = 0;                            % liquid volume at previous temp. step
    prevV_l(2:nrows)    = COMP.V_l{i}(1:nrows-1);
    prevrho_l           = [];
    prevrho_l(1)        = 0;                            % liquid density at previous temp. step
    prevrho_l(2:nrows)  = COMP.rho_l{i}(1:nrows-1);
    
    dT                  = [];
    dT(1)               = 0;                            % change in temperature
    dT(2:nrows)         = COMP.temps{i}(2:nrows) - prevT(2:nrows);   
    dV_l                = [];
    dV_l(1)             = 0;                            % change in liquid volume
    dV_l(2:nrows)       = COMP.V_l{i}(2:nrows) - prevV_l(2:nrows);
    drho_l              = [];
    drho_l(1)           = 0;                            % change in liquid density
    drho_l(2:nrows)     = COMP.rho_l{i}(2:nrows) - prevrho_l(2:nrows);
    
    % energy variation in liquid:
    dE_cool             = [];
    dE_cool(1)          = 0;                           % energy spent in cooling    
    dE_cool(2:nrows)    = dT(2:nrows).*prevrho_l(2:nrows)*RES.cp.*prevV_l(2:nrows);
    dE_freez            = [];
    dE_freez(1)         = 0;                           % energy spent in freezing    
    dE_freez(2:nrows)   = dV_l(2:nrows).*prevrho_l(2:nrows)*ICE.L;                      
    dE                  = [];
    dE                  = dE_cool + dE_freez;          % total d energy
    
    % cumulative d energies
    DeltaE_cool         = [];
    DeltaE_cool         = cumsum(dE_cool);
    DeltaE_freez        = [];
    DeltaE_freez        = cumsum(dE_freez);
    DeltaE              = [];
    DeltaE              = cumsum(dE);
    COMP.DeltaE{i}      = DeltaE;
    COMP.fE_rmn{i}      = (E_init+DeltaE)/E_init;       % fraction of energy remaining

% end

disp('Composition data imported')

end
