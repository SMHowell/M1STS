%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function COMP = getSalts(COMP, i, nrows)

try
    COMP.Ca{i}     = table2array(COMP.Table{i}(:,["Ca (mol/(kg water))"])).';
catch
    COMP.Ca{i}     = zeros(nrows,1);
end

try
    COMP.Mg{i}     = table2array(COMP.Table{i}(:,["Mg (mol/(kg water))"])).';
catch
    COMP.Mg{i}     = zeros(nrows,1);
end

try
    COMP.Na{i}     = table2array(COMP.Table{i}(:,["Na (mol/(kg water))"])).';
catch
    COMP.Na{i}     = zeros(nrows,1);
end    
    
try    
    COMP.K{i}      = table2array(COMP.Table{i}(:,["K (mol/(kg water))"])).';
catch
    COMP.K{i}      = zeros(nrows,1);
end   
    
try    
    COMP.Cl{i}     = table2array(COMP.Table{i}(:,["Cl (mol/(kg water))"])).';
catch
    COMP.Cl{i}     = zeros(nrows,1);
end
    
try    
    COMP.S{i}      = table2array(COMP.Table{i}(:,["S (mol/(kg water))"])).';
catch
    COMP.S{i}      = zeros(nrows,1);
end

try
    COMP.C{i}      = table2array(COMP.Table{i}(:,["C (CO2+carbonate) (mol/(kg water))"])).';
catch
    COMP.C{i}      = zeros(nrows,1);
end

try
    COMP.Si{i}     = table2array(COMP.Table{i}(:,["Si (mol/(kg water))"])).';
catch
    COMP.Si{i}     = zeros(nrows,1);
end

try
    COMP.Mtg{i}    = table2array(COMP.Table{i}(:,["Mtg (mol/(kg water))"])).';
catch
    COMP.Mtg{i}    = zeros(nrows,1);
end


end