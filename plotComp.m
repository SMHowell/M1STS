%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)


semilogy(OUT.t2(1:IN.outInd2)/365.25/24/3600,[OUT.comp{1,:}],'LineWidth',2);
hold on;
semilogy(OUT.t2(1:IN.outInd2)/365.25/24/3600,[OUT.comp{2,:}],'LineWidth',2);
hold on;
semilogy(OUT.t2(1:IN.outInd2)/365.25/24/3600,[OUT.comp{3,:}],'LineWidth',2);
hold on;
semilogy(OUT.t2(1:IN.outInd2)/365.25/24/3600,[OUT.comp{4,:}],'LineWidth',2);
hold on;
semilogy(OUT.t2(1:IN.outInd2)/365.25/24/3600,[OUT.comp{5,:}],'LineWidth',2);
hold on;
semilogy(OUT.t2(1:IN.outInd2)/365.25/24/3600,[OUT.comp{6,:}],'LineWidth',2);
hold on;
semilogy(OUT.t2(1:IN.outInd2)/365.25/24/3600,[OUT.comp{4,:}],'LineWidth',2);
hold on;
semilogy(OUT.t2(1:IN.outInd2)/365.25/24/3600,[OUT.comp{8,:}],'LineWidth',2);

grid on;
% for i = 1:IN.nErupt
    % xline(OUT.eruptTimes(i))
% end

legend('Ca', 'Mg', 'Na', 'K', 'Cl', 'S', 'C', 'Si');