%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function COMP = initializeEnergy(IN,COMP,M)

% Calculate the energy loss from liquid at 
% each temperature step and how it's spent

i = IN.simu;
nrows = length(COMP.temps{i});

% total energy in reservoir at liquidus
E_init              = COMP.rho_l{i}(1)*COMP.V_init{i}(1) * (COMP.temps{i}(1)*M.CpRes + IN.L);

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
dE_cool(2:nrows)    = dT(2:nrows).*prevrho_l(2:nrows)*M.CpRes.*prevV_l(2:nrows);
dE_freez            = [];
dE_freez(1)         = 0;                           % energy spent in freezing    
dE_freez(2:nrows)   = dV_l(2:nrows).*prevrho_l(2:nrows)*IN.L;                      
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

end