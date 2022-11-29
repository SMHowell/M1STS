%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MISC] = polarPlotInt(IN,M,OUT,MISC)

%%%%%%%%%%%%%%%%%%%%%%%
% Generate plot if engaged
%%%%%%%%%%%%%%%%%%%%%%%
if IN.pltOn
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Initialize plot, axes, movie frames
    %%%%%%%%%%%%%%%%%%%%%%%
    if M.step==1
        figNum = 10;
        
        % Set up for movie if checked
        if IN.movOn
            % Plotting needs to be on
            IN.pltOn = 1;
            
            % Set up counters and arrays
            MISC.F1  = struct('cdata', cell(1,IN.outN), 'colormap', cell(1,IN.outN));
            MISC.fg  = figure(figNum);
            set(figure(figNum),'Units','pixels','Position',[100 100 1920 1080])
            set(figure(figNum),'doublebuffer','off','Visible','Off');
        else
            MISC.fg = figure(figNum);
            set(figure(figNum),'Units','pixels','Position',[100 100 1920 1080])
        end
    end
    
    % Select current figure without stealing focus
    set(0,'CurrentFigure',MISC.fg); clf; set(gcf,'color','k');
    fontSize = 24;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Polar plot
    %%%%%%%%%%%%%%%%%%%%%%%
    % Get reservoir geometry
    rPol = M.r/1e3;
    
    %%%%%%%%%%%%%
    % Cut out reservoir from temperature structure for plotting away from
    % reservoir
    % Polar arrays
    sp2    = subplot(122);
    Tpol   = M.T;  % Radial Temperatue
    width  = pi;   % Radial span (radians)
    theta  = linspace(pi - width/2,pi + width/2,length(rPol)); % Polar angle

    % Create plottable grids from 1D arrays
    [~,R2] = meshgrid(theta,rPol);      % Radius
    T2     = interp1(rPol,Tpol,R2);     % Temperature
    Xpol   = rPol'*cos(theta);          % X-coord on plot
    Ypol   = rPol'*sin(theta); % Y-coord on plot

    % Main plot (no reservoir thermal effects
    hold off;
    h = pcolor(Xpol,Ypol,T2); 

    xlabel('Radial Distance [km]');
    ylabel('Radial Distance [km]');
    set(h,'linestyle','none'); axis equal; shading interp;
    zRnd = 1e3*ceil(M.z(1)/1e6);
    axis([-zRnd, zRnd, -zRnd, zRnd]);
    
    %%% Boundaries
    hold on; alpha = 0.85;
    
    % Ice shell
    color = 'cyan';
    r1 = M.r(end-1)/1e3; 
    r2 = M.rOcn/1e3;
    xBnd1 = r1 * cos(theta+pi); yBnd1 = r1 * sin(theta+pi); 
    xBnd2 = r2 * cos(theta+pi); yBnd2 = r2 * sin(theta+pi); 
    fill([xBnd1,xBnd2], [yBnd1,yBnd2],color,'facealpha',alpha);
    plot(-xBnd1, -yBnd1,'color','k');
    
    % Ocean
    color = [0.5843    0.8157    0.9882];
    r1 = M.rOcn/1e3; 
    r2 = M.rSil/1e3;
    xBnd1 = r1 * cos(theta+pi); yBnd1 = r1 * sin(theta+pi); 
    xBnd2 = r2 * cos(theta+pi); yBnd2 = r2 * sin(theta+pi); 
    fill([xBnd1,xBnd2], [yBnd1,yBnd2],color,'facealpha',alpha);
    plot(-xBnd1, -yBnd1,'color','k');
    
    % Mantle
    color = [0.7961    0.3765    0.0824];
    r1 = M.rSil/1e3; 
    r2 = M.rIrn/1e3;
    xBnd1 = r1 * cos(theta+pi); yBnd1 = r1 * sin(theta+pi); 
    xBnd2 = r2 * cos(theta+pi); yBnd2 = r2 * sin(theta+pi); 
    fill([xBnd1,xBnd2], [yBnd1,yBnd2],color,'facealpha',alpha);
    plot(-xBnd1, -yBnd1,'color','k');
    
    % Core
    color = [1 1 1]/4;
    r1 = M.rIrn/1e3; 
    xBnd1 = r1 * cos(theta+pi); yBnd1 = r1 * sin(theta+pi); 
    fill(xBnd1, yBnd1,color,'facealpha',alpha);
    plot(-xBnd1, -yBnd1,'color','k');
    
    %%% Make it Pretty
    set(gca,'fontsize',fontSize)
    colormap(parula(256));

    cb = colorbar;
    ylabel(cb,'Temperature [K]','fontsize',fontSize);
    set(cb,'xcolor','w','ycolor','w');
    caxis([100,1300]);
    
    set(gca,'layer','top','color','k','xcolor','w','ycolor','w');
    title([num2str(MISC.outTime ,'%05.1f') ' '  MISC.unit],'fontname','Monowidth','color','w');
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Temperature 
    %%%%%%%%%%%%%%%%%%%%%%%
    sp1 = subplot(121);
    Tmax = 2000;
        
    %%% Boundaries
    hold on; alpha = 0.35;
    
    % Ice
    color = 'cyan';
    r1 = M.r(end-1)/1e3; r2 = M.rOcn/1e3;     
    fill([0,Tmax,Tmax,0,0], [r2,r2,r1,r1,r2],color,'facealpha',alpha)
    fill([0,Tmax,Tmax,0,0],-[r2,r2,r1,r1,r2],color,'facealpha',alpha)
    
    % Ocean
    color = [0.5843    0.8157    0.9882];
    r1 = M.rOcn/1e3; r2 = M.rSil/1e3;     
    fill([0,Tmax,Tmax,0,0], [r2,r2,r1,r1,r2],color,'facealpha',alpha)
    fill([0,Tmax,Tmax,0,0],-[r2,r2,r1,r1,r2],color,'facealpha',alpha) 
    
    % Mantle
    color = [0.7961    0.3765    0.0824];
    r1 = M.rSil/1e3; r2 = M.rIrn/1e3;     
    fill([0,Tmax,Tmax,0,0], [r2,r2,r1,r1,r2],color,'facealpha',alpha)
    fill([0,Tmax,Tmax,0,0],-[r2,r2,r1,r1,r2],color,'facealpha',alpha) 
    
    % Core
    color = 'white';
    r1 = M.rIrn/1e3; r2 = -r1;     
    fill([0,Tmax,Tmax,0,0], [r2,r2,r1,r1,r2],color,'facealpha',alpha)
    
    %%% Temperature
    plot([fliplr(Tpol),Tpol],[-fliplr(rPol),rPol],'w','linewidth',2);
    
    %%% Make it pretty
    xlabel('Temperature [K]');
    ylabel('Depth [km]');
    axis([0, Tmax,-zRnd,zRnd]);
    set(gca,'layer','top','color','k','xcolor','w','ycolor','w');
    set(gca,'fontsize',fontSize); grid on; box on;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Better Align
    %%%%%%%%%%%%%%%%%%%%%%%
    
    set(sp1,'position',[0.10,0.15,0.233,0.75]);
    set(sp2,'position',[0.45,0.15,2*0.233,0.75]);
    drawnow;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Movie Frames
    %%%%%%%%%%%%%%%%%%%%%%%
    
    if IN.movOn
        MISC.F1(IN.outInd) = getframe(gcf);
    end
    
    % Save movie to file
    if (IN.outInd == IN.NOut) && IN.movOn
        noFrame = zeros(IN.outInd,1);
        for i=1:IN.outInd
            if isempty(MISC.F1(i).cdata)
                noFrame(i) = 1;
            end
        end
        MISC.F1(logical(noFrame)) = [];
        if IN.resOn
        outputName = ['H0' num2str(IN.H0_ice/1e3,'%0.1f') '_rRes' num2str(IN.rRes/1e3,'%0.1f') '_zRes' num2str(IN.zResTop/1e3,'%0.1f')];
        else
               outputName = ['H0' num2str(IN.H0_ice/1e3,'%0.1f') '_resOff'];    
        end
        v = VideoWriter(['Output\Movies\' outputName '.mp4'],'MPEG-4');
        if IN.outInd >= 300
            v.FrameRate = 60;
        elseif IN.outInd >= 150
            v.FrameRate = 30;
        else
            v.FrameRate = 12;
        end
        v.Quality   = 100;
        open(v); writeVideo(v,MISC.F1); close(v)
        clear MISC;
        MISC = [];
    end % End movie save
end % end plot
end % End function








