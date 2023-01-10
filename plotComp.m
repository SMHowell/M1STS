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
% hold on; 

% Get last index populated before too many eruptions occurred
lastIndTime = find(OUT.t2>0,1,'last');

% Colormap
cc = parula(8);

for i=1:IN.Ncomp
    semilogy((OUT.t2(1:lastIndTime)-IN.tRes)/IN.kyr2s,OUT.comp(i,1:lastIndTime), 'color', cc(i,:), 'linewidth',2);
    hold on;
end

for i=1:50:M.eruption-1
    if OUT.eruptTimes(i) < 110 * IN.kyr2s
        xline([(OUT.eruptTimes(i)-IN.tRes)/IN.kyr2s,(OUT.eruptTimes(i)-IN.tRes)/IN.kyr2s],'k', 'linewidth', 1.5);
        hold on;
    end
end

legend(M.compLabels,'location','eastoutside');

axis tight; box on; grid on; 
xlim([0,OUT.t2(lastIndTime)/IN.kyr2s])
xlabel('Time [kyr]')
ylabel('Concentration [mol/kg H_2O]')
set(gca,'fontsize',20);





figure(3); clf; set(gcf,'color','w'); 

%  eruption intervals and volumes
subplot(221);

yyaxis left
plot(OUT.eruptTimes(1:M.eruption-2)/IN.kyr2s,diff(OUT.eruptTimes(1:M.eruption-1)/IN.yr2s), 'LineWidth', 1.5)

xlim([0,OUT.t2(lastIndTime)/IN.kyr2s])
xlabel('Time [kyr]')
ylabel('Eruption Interval [yr]')

yyaxis right

% xPlot = OUT.t2(OUT.eruptV~=0)/IN.kyr2s;
% yPlot = OUT.eruptV(OUT.eruptV~=0)/1e9;
xPlot = OUT.t2(1:length(OUT.eruptV))/IN.kyr2s;
yPlot = cumsum(OUT.eruptV)/1e9;
ylim([0,0.12])

plot(xPlot(2:end),yPlot(2:end), 'LineWidth', 1.5)
xlim([0,OUT.t2(lastIndTime)/IN.kyr2s])
xlabel('Time [kyr]')
ylabel('Eruption Volume [km^3]')

box on; grid on; 
set(gca,'fontsize',20);

% Tmelt
subplot(222);

xPlot = OUT.t(OUT.Tmelt~=0)/IN.kyr2s;
yPlot = OUT.Tmelt(OUT.Tmelt~=0);

plot(xPlot(2:end),yPlot(2:end), 'k', 'LineWidth', 1.5)

box on; grid on; 
xlim([0,OUT.t2(lastIndTime)/IN.kyr2s])
xlabel('Time [kyr]')
ylabel('Reservoir Freezing Temperature [K]')
set(gca,'fontsize',20);


% Ice thickness 
subplot(223);

xPlot = OUT.t(OUT.Dice~=0)/IN.kyr2s;
yPlot = OUT.Dice(OUT.Dice~=0)/1e3;

plot(xPlot(2:end),yPlot(2:end), 'k', 'LineWidth', 1.5)

box on; grid on; 
xlim([0,OUT.t2(lastIndTime)/IN.kyr2s])
xlabel('Time [kyr]')
ylabel('Local Ice Shell Thickness [km]')
set(gca,'fontsize',20);


% Surface Temp
subplot(224);

xPlot = OUT.t(OUT.Dice~=0)/IN.kyr2s;
yPlot = OUT.Tsurf(OUT.Tsurf~=0)-OUT.Tsurf(1);

plot(xPlot(2:end),yPlot(2:end), 'k', 'LineWidth', 1.5)

box on; grid on; 
xlim([0,OUT.t2(lastIndTime)/IN.kyr2s])
xlabel('Time [kyr]')
ylabel('Surface Temperature [K]')
set(gca,'fontsize',20);

