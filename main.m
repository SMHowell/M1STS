%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%******** RUN FROM resInput.m **********

%%%%%%%%%%%%%%%%%%
% Initialize Model
%%%%%%%%%%%%%%%%%%
initialize;

%%%%%%%%%%%%%%%%%%
% Main Time Loop
%%%%%%%%%%%%%%%%%%
while M.t < IN.tMax
    
    % Get timestep
    M.dt = getTimestep(M);
        
    % Temperature diffusion solve
    M = thermalSolver(M,IN,BOD);
        
    % Move interfaces
    M = meltingFreezing(M,IN);
    
    % Update time 
    M.t    = M.t + M.dt;
    M.step = M.step + 1;
    
    % Determine whether to emplace reservoir
    M = reservoirEmplacement(M,IN);
    
    % Update thermal properties
    M = getThermalProperties(M);
    
    % Generate Output
    [IN, OUT, MISC] = outputManager(M,IN,BOD,OUT,MISC);

       
end