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
    set(0,'CurrentFigure',MISC.fg); clf; whitebg('k'); set(gcf,'color','k');
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
    theta  = linspace(- width/2,width/2,length(rPol)); % Polar angle

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
    set(h,'linestyle','none'); axis equal tight; shading interp;
        
    %%% Boundaries
    hold on;
    
    % Surface
    r    = M.r(end-1)/1e3;
    xBnd = r * cos(theta); yBnd = r * sin(theta); h = plot(xBnd, yBnd);
    set(h,'color','w','linewidth',2)
    
    % Ice shell
    r    = M.rOcn/1e3;
    xBnd = r * cos(theta); yBnd = r * sin(theta); h = plot(xBnd, yBnd);
    set(h,'color','w','linewidth',2)
    
    % Mantle
    r    = M.rSil/1e3;
    xBnd = r * cos(theta); yBnd = r * sin(theta); h = plot(xBnd, yBnd);
    set(h,'color','w','linewidth',2)
    
    % Core
    r    = M.rIrn/1e3;
    xBnd = r * cos(theta); yBnd = r * sin(theta); h = plot(xBnd, yBnd);
    set(h,'color','w','linewidth',2)
       
    %%% Make it Pretty
    set(gca,'fontsize',fontSize)
    colormap(parula(256));

    cb = colorbar;
    ylabel(cb,'Temperature [K]','fontsize',fontSize)
    
    set(gca,'layer','top')
    title([num2str(MISC.outTime ,'%05.1f') ' '  MISC.unit],'fontname','Monowidth');
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Temperature and Frozen Fraction
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Temperature
    sp1 = subplot(121);
    plot([fliplr(Tpol),Tpol],[-fliplr(rPol),rPol]);
    xlabel('Temperature [K]');
    ylabel('Depth [km]');
    axis([0, 3000,-M.z(1)/1e3,M.z(1)/1e3]);
    set(gca,'fontsize',fontSize); grid on;box on;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Better Align
    %%%%%%%%%%%%%%%%%%%%%%%
    
    set(sp1,'position',[0.08,0.10+(0.8-0.55)/2,0.233,0.55]);
    set(sp2,'position',[0.38,0.10+(0.8-0.55)/2,0.233,0.55]);
    set(sp3,'position',[0.73,0.10+(0.8-0.55)/2,0.233,0.55]);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Movie Frames
    %%%%%%%%%%%%%%%%%%%%%%%
    
    if IN.movOn
        MISC.F1(IN.outInd) = getframe(gcf);
    end
    drawnow
    
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








