%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()

hold on
plot(OUT.t2(1:IN.outInd2),[OUT.comp{1,:}])
plot(OUT.t2(1:IN.outInd2),[OUT.comp{2,:}])
plot(OUT.t2(1:IN.outInd2),[OUT.comp{3,:}])
plot(OUT.t2(1:IN.outInd2),[OUT.comp{4,:}])
plot(OUT.t2(1:IN.outInd2),[OUT.comp{5,:}])
plot(OUT.t2(1:IN.outInd2),[OUT.comp{6,:}])
plot(OUT.t2(1:IN.outInd2),[OUT.comp{4,:}])
plot(OUT.t2(1:IN.outInd2),[OUT.comp{8,:}])

i=0
for i in IN.nErupt
    
    i = i+1
end


hold off

legend('Ca', 'Mg', 'Na', 'K', 'Cl', 'S', 'C', 'Si');