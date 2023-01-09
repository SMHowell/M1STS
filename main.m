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
rOld =  BOD.R-M.rOcn;
drmax = 0;
while M.t < IN.tMax
    rNew = BOD.R-M.rOcn;
    if abs(rNew-rOld)>drmax
        drmax = abs(rNew-rOld);
        disp(['Step ' num2str(M.step) ', dr = ' num2str(drmax)]);
    end
%     if M.step == 58413
%         fasdf
%     end
    rOld = rNew;
    
    
    % Generate Output
    [IN, OUT, MISC] = outputManager(M,IN,BOD,OUT,MISC);
    
    % Get timestep
    M.dt = getTimestep(M);
    
    % Temperature diffusion solve
    M = thermalSolver(M,BOD,IN,MAT);
    
    % Differentiation
    [M,OUT] = differentiate(M,IN,BOD,MAT,OUT);
    
    % Move interfaces
    M = meltingFreezing(M,IN,BOD,MAT);
    
    % Update time 
    M.t    = M.t + M.dt;
    M.step = M.step + 1;
    
    % Update thermal properties
    M = getThermalProperties(M,IN,MAT);

    % Check for convection
    M = getConvection(M,MAT);
       
end