%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up fig
figure(2); clf; whitebg('w'); set(gcf,'color','w');  
hold on; 

% Get last index populated before too many eruptions occurred

lastIndTime = find(OUT.t2>0,1,'last');

for i=1:IN.Ncomp
    plot(OUT.t2(1:lastIndTime)/IN.kyr2s,OUT.comp(i,1:lastIndTime),'linewidth',2);
end
legend(M.compLabels,'location','eastoutside');

for i=1:25:M.eruption-1
    plot([OUT.eruptTimes(i)/IN.kyr2s,OUT.eruptTimes(i)/IN.kyr2s],[0,1],'r');
end

axis tight; box on; grid on; 
xlim([0,OUT.t2(lastIndTime)/IN.kyr2s])
xlabel('Time [kyr]')
ylabel('Concentration [mol/kg H_2O]')
set(gca,'fontsize',24);


% Now eruption intervals
figure(3); clf; set(gcf,'color','w');  

plot(OUT.eruptTimes(1:M.eruption-2)/IN.kyr2s,diff(OUT.eruptTimes(1:M.eruption-1)/IN.yr2s))

axis tight; box on; grid on; 
xlim([0,OUT.t2(lastIndTime)/IN.kyr2s])
xlabel('Time [kyr]')
ylabel('Eruption Interval [yr]')
set(gca,'fontsize',24);