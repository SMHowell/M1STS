%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MISC] = polarPlot(IN,M,OUT,MISC)

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
    yRef  = M.r(end-1);               % Surface reference depth to Y-scale plot
    rPol  = M.r(M.iOcnTop:end-1)/1e3; % Polar plot radii
    rResCenter_0 = (IN.zResTop+IN.rRes)/1e3; % initial center of reservoir
    rRes_0 = IN.rRes/1e3;             % initial radius of reservoir
    
    rResCenter = mean(yRef-[M.rResTop,M.rResBot])/1e3;
    rRes       = M.rRes/1e3;
    
    %%%%%%%%%%%%%
    % Cut out reservoir from temperature structure for plotting away from
    % reservoir
    % Polar arrays
    sp2 = subplot(132);
    Tpol  = M.T(M.iOcnTop:end-1);     % Radial Temperatue
    Tpol  = Tpol(end) * (Tpol(1)/Tpol(end)).^((rPol(end)-rPol)/(rPol(end)-rPol(1)));  % Masked reservoir w/ stand in conductive temperature structure
    theta  = linspace(pi/2-0.05,pi/2+0.05,length(rPol)); % Polar angle
    
    % Create plottable grids from 1D arrays
    [~,R2] = meshgrid(theta,rPol);      % Radius
    T2     = interp1(rPol,Tpol,R2);     % Temperature
    Xpol   = rPol'*cos(theta);          % X-coord on plot
    Ypol   = yRef/1e3-rPol'*sin(theta); % Y-coord on plot
    
    % Main plot (no reservoir thermal effects
    h = pcolor(Xpol,Ypol,T2); hold on;
    
    %%%%%%%%%%%%%
    % Plot local reservoir effects
    % Polar arrays
    if M.vRes>0
        Tpol  = M.T(M.iOcnTop:end-1);     % Radial Temperatue
        resTheta = M.rRes/(M.r(end-1)-rResCenter*1e3); % Angle subtended by reservoir
        theta = linspace(pi/2-resTheta,pi/2+resTheta,length(rPol)); % Polar angle
        
        % Create plottable grids from 1D arrays
        [~,R2] = meshgrid(theta,rPol);      % Radius
        T2     = interp1(rPol,Tpol,R2);     % Temperature
        Xpol   = rPol'*cos(theta);          % X-coord on plot
        Ypol   = yRef/1e3-rPol'*sin(theta); % Y-coord on plot
        
        h = pcolor(Xpol,Ypol,T2);
    end
    xlabel('Depth [km]');
    ylabel('Distance [km]');
    set(h,'linestyle','none'); axis equal tight; shading interp;
    set(gca,'ydir','reverse','fontsize',fontSize)
    colormap(parula(256));
    
    axis([-15, 15, -5, 35]);
    
    cb = colorbar;
    ylabel(cb,'Temperature [K]','fontsize',fontSize)
    
    title([num2str(MISC.outTime ,'%05.1f') ' '  MISC.unit],'fontname','Monowidth');
    
    %%%%%%%%%%%%%
    % Plot the reservoir
    % Plot outline of original reservoir and outline of current config
    p = nsidedpoly(100, 'Center', [0 rResCenter_0], 'Radius', rRes_0);
    plot(p, 'FaceColor', 'c','edgecolor','c')
    
    if M.rRes>0
        p = nsidedpoly(100, 'Center', [0 rResCenter], 'Radius', rRes);
        plot(p, 'FaceColor', 'c','edgecolor','k')
        set(gca,'layer','top');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Temperature and Frozen Fraction
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % Temperature
    sp1 = subplot(131);
    plot(M.T,M.z/1e3);
    set(gca,'ydir','reverse');
    xlabel('Temperature [K]');
    ylabel('Depth [km]');
    axis([95, 300,0,30]);
    set(gca,'fontsize',fontSize); grid on;box on;
    
    % Frozen Fraction
    sp3 = subplot(133);
    plot(OUT.t(1:IN.outInd)/IN.Myr2s,OUT.fV(1:IN.outInd));
    xlabel('Time [Myr]');
    ylabel('Frozen Volume Fraction');
    axis([0,IN.tMax/IN.Myr2s,0,1]);
    set(gca,'fontsize',fontSize); grid on; box on;
    
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
        outputName = ['H0' num2str(IN.H0/1e3,'%0.1f') '_rRes' num2str(IN.rRes/1e3,'%0.1f') '_zRes' num2str(IN.zResTop/1e3,'%0.1f')];
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








