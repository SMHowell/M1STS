%% Composite Iron Rheology
clear; clc; 

% Input params
figures = 1; % 1 = Turn on figures
T_lo = 900; % Lower Temp Bound [K]
T_m = 1810; % Melting Temp / Upper Bound [K]

T = linspace(T_lo , T_m, 500); % Temperatures of interest [K]
R = 8.314; % Ideal gas constant [kJ mol^-1 K^-1]
d = logspace(-7,1,500); % Grain size [m]

w_I = 4.11e-5; % Io orbital forcing frequency [1/s]
w_E = 2.05e-5; % Europa orbital forcing frequency [1/s]
w_G = 1.02e-5; % Ganymede orbital forcing frequency [1/s]

% Make arrays
[T2,d2] = meshgrid(T,d);
% load VanceOut;
% V_imk2 = k2_output(:,2);
% V_eta = 10.^(k2_output(:,1));

%% Rheological properties - alpha iron
% Diffusion creep
Vmol    = 7.11e-6;  % Molar volume [m^3]
Dov     = 2.0e-4;   % Volume diffusion coeff. [m^2/s]
Dob     = 2.2e-3;   % Grain boundary diffusion coeff. [m^2/s]
Eav     = 251e3;    % Volume activation energy [J/mol]
Eab     = 174e3;    % GB activation energy [J/mol]
delta   = 4.96e-10; % GB width [m]

% Prefactors
Dv = Dov .* exp( - Eav ./ (R .*(T2)));
Db = Dob .* exp( - Eab ./ (R .*(T2)));

% Full Newtonian rheology form from Goldsby and Kohlstedt 2001
eta_a = ((42 .* Vmol ./ (R .*(T2) .* ( d2.^2))).*...
    (Dv + pi *  delta .*   Db ./ d2)).^-1;

%% Rheological properties - gamma iron
% Diffusion creep
delta   = 5.16e-10; % GB width [m]
Vmol    = 7.29e-6;  % Molar volume [m^3]
Dov     = 1.8e-4;   % Volume diffusion coeff. [m^2/s]
Dob     = 7.5e-14/delta; % Grain boundary diffusion coeff. [m^2/s]
Eav     = 270e3;    % Volume activation energy [J/mol]
Eab     = 159e3;    % GB activation energy [J/mol]

% Prefactors
Dv = Dov .* exp( - Eav ./ (R .*(T2)));
Db = Dob .* exp( - Eab ./ (R .*(T2)));

% Full Newtonian rheology form from Goldsby and Kohlstedt 2001
eta_g = ((42 .* Vmol ./ (R .*(T2) .* ( d2.^2))).*...
    (Dv + pi *  delta .*   Db ./ d2)).^-1;

%% Rheological properties - delta iron
% Diffusion creep
% Diffusion creep
Vmol    = 7.11e-6;  % Molar volume [m^3]
Dov     = 1.9e-4;   % Volume diffusion coeff. [m^2/s]
Dob     = 1.1e-12/delta; % Grain boundary diffusion coeff. [m^2/s]
Eav     = 239e3;    % Volume activation energy [J/mol]
Eab     = 174e3;    % GB activation energy [J/mol]
delta   = 4.96e-10; % GB width [m]

% Prefactors
Dv = Dov .* exp( - Eav ./ (R .*(T2)));
Db = Dob .* exp( - Eab ./ (R .*(T2)));

% Full Newtonian rheology form from Goldsby and Kohlstedt 2001
eta_d = ((42 .* Vmol ./ (R .*(T2) .* ( d2.^2))).*...
    (Dv + pi *  delta .*   Db ./ d2)).^-1;

%% Now create full rheology
T_a = 1184; % Max alpha iron temp [K]
T_g = 1665; % Max gamma iron temp [K]

eta = zeros(size(T2)); % Full viscosity [Pa s]
eta(T2<=T_a) = eta_a(T2<=T_a);
eta(T2>T_a & T2<=T_g) = eta_g(T2>T_a & T2<=T_g);
eta(T2>T_g) = eta_d(T2>T_g);


%% Shear Modulus
T_ref = 300; % Reference Shear modulus temperature [K]

% Refence modulus
G0_a = 64e9; % alpha iron
G0_g = 81e9; % gamma iron
G0_d = 39e9; % delta iron

% Tm/G0 * dGmod/dT
depT_a = -0.81; % alpha iron
depT_g = -0.91; % gamma iron
depT_d = -0.72; % delta iron

% Gmod(T)
G_a = G0_a * (1 + ((T2-T_ref)/T_m) *depT_a); % alpha iron
G_g = G0_g * (1 + ((T2-T_ref)/T_m) *depT_g); % gamma iron
G_d = G0_d * (1 + ((T2-T_ref)/T_m) *depT_d); % delta iron

Gmod = zeros(size(T2)); % Full Modulus [Pa]
Gmod(T2<=T_a) = G_a(T2<=T_a);
Gmod(T2>T_a & T2<=T_g) = G_g(T2>T_a & T2<=T_g);
Gmod(T2>T_g) = G_d(T2>T_g);

% Dissipation
q_e2_I = w_I^2 * eta ./ (1 + w_I^2 * eta.^2./Gmod.^2);
q_e2_E = w_E^2 * eta ./ (1 + w_E^2 * eta.^2./Gmod.^2);
q_e2_G = w_G^2 * eta ./ (1 + w_G^2 * eta.^2./Gmod.^2);

%% Find Major Contributors to Global Heat flux
didx = 1e-3; % Grain size of interest
[v,idx] = min(abs(d-didx)); % index of grain size

eta_d = eta(idx,:); % Viscosity at one grain size
G_d   = Gmod(idx,:); % Modulus at one grain size

%% Silicate Viscosity
% Hirth and Kohlstedt Wet Peridotite
A_Si_H2O = 1e24;   % Viscosity Prefactor, hydrated [Pa s * (Pa/m^3)]

p_Si = 3;      % Stress exponent
fH2O = 1e9;    % Water Fugacity [Pa] -- may be very large for hydrated silicates  
               % https://people.earth.yale.edu/sites/default/files/files/Karato/32Otsuka-Karato%20(2011).pdf
E_Si = 375e3;  % Activation energy [J/kg K]

eta_Si_H2O   = A_Si_H2O * (d2.^(p_Si)/fH2O) .* exp(E_Si ./ (R * T2)); % Viscosity
eta_Si_H2O_d = eta_Si_H2O(idx,:); % Viscosity at one grain size

%% Figure
figure(1); clf; set(gcf,'color','w');

% Wet
h = pcolor(T2,d2,eta./eta_Si_H2O); cb = colorbar;
set(h,'linestyle','none'); shading interp;

xmin = min(T2(:)); xmax = max(T2(:));
ymin = 1e-5; ymax = 1;   ylim([ymin,ymax]); 
cmin = 1e-12; cmax = 1e-2; caxis([cmin, cmax]);

set(gca,'ColorScale','log','yscale','log');
set(gca,'ytick',[1e-5,1e-4,1e-3,1e-2,1e-1,1])
set(gca,'yticklabel',{'10 \Gmodm','100 \Gmodm','1 mm','1 cm','10 cm','1 m'})

ylabel(cb,'\eta_{Fe}/ \eta_{Si}')
xlabel('Core Temperature [K]');
ylabel('Grain Size');

grid on; set(gca, 'layer', 'top'); set(gca,'fontsize',24); axis square;

%% Figure
figure(2); clf; set(gcf,'color','w');

fFe       = linspace(0,1,500); % Volume fraction iron
[f2, Tc2] = meshgrid(fFe,T);

eta2_Si   = interp1(T,eta_Si_H2O_d,Tc2);
eta2_Fe   = interp1(T,eta_d,Tc2);

eta_comp = (f2./eta2_Fe + (1-f2)./eta2_Si).^-1; % Composite viscosity
h = pcolor(Tc2,f2,eta_comp);
set(h,'linestyle','none'); 

cb = colorbar;
ylabel(cb,'Composite Viscosity [Pa s]')

xlabel('Temperature [K]')
ylabel('Iron Volume Fraction')


axis square; shading interp;
set(gca,'ColorScale','log');
caxis([1e14,1e20]);
box on; grid on; set(gca,'fontsize',24); set(gca, 'layer', 'top'); 

%% Plot figure for paperi had a momentum an
Period = 2*pi*eta./Gmod; % Orbital period that maximizes dissipation [s]

% Set orbital periods of moons
P_I = 1.769;
P_E = 3.551;
P_G = 7.115;
P_C = 16.69;
P_T = 5.877;

% Make figure
figure(3); clf; set(gcf,'color','w');

h = pcolor(T2,d2,Period/(3600*24)); cb = colorbar;
set(h,'linestyle','none'); shading interp;

xmin = min(T2(:)); xmax = max(T2(:));
ymin = 1e-5; ymax = 1;   ylim([ymin,ymax]); 
cmin = 1/24; cmax = 1e4; caxis([cmin, cmax]);
% Slowest is Nesso at Neptune (27 years)

cb.Ticks = [1/24, 1, 30.4, 365e0, 365e1, 1e4];
cb.TickLabels = {'1 hour','1 day','1 month','1 year','10 years','26.7 years'};

set(gca,'ColorScale','log','yscale','log');
set(gca,'ytick',[1e-5,1e-4,1e-3,1e-2,1e-1,1])
set(gca,'yticklabel',{'10 \Gmodm','100 \Gmodm','1 mm','1 cm','10 cm','1 m'})

xlabel('Core Temperature [K]');
ylabel('Core Grain Size');

grid on; set(gca, 'layer', 'top'); set(gca,'fontsize',24); axis square;

% Remove weird periods
hold on; 

[c1, ~] = contour(T2,d2,Period/(3600*24),[cmin, cmin],'w','linewidth',1);
[c2, ~] = contour(T2,d2,Period/(3600*24),[cmax, cmax],'w','linewidth',1);

x1 = c1(1,2:end);
y1 = c1(2,2:end);

x2 = c2(1,2:end);
y2 = c2(2,2:end);

fill([x1, xmax, xmin],[y1, ymin, ymin],'w','edgecolor','w')
fill([x2, xmax, xmin],[y2, ymax, ymax],'w','edgecolor','w')

% Add body labels

contour(T2,d2,Period/(3600*24),[P_I, P_I],'k','linewidth',2);
contour(T2,d2,Period/(3600*24),[P_E, P_E],'r','linewidth',2);
contour(T2,d2,Period/(3600*24),[P_G, P_G],'k','linewidth',2);
contour(T2,d2,Period/(3600*24),[P_C, P_C],'k','linewidth',2);
contour(T2,d2,Period/(3600*24),[P_T, P_T],'k','linewidth',2);



% %% Plots
% if figures
% %%%%%%%%%%%%%%%%%
% % Viscosity map
% %%%%%%%%%%%%%%%%%
% figure(1); clf; set(gcf,'color','w'); fontsize = 16;
% labels = {'0.1 um','1 um','10 um','0.1 mm','1 mm','1 cm','10 cm','1 m','10 m'};
% 
% subplot(131)
% [C,h] = contourf(T2,log10(d2),log10(eta),[1:25]);
% clabel(C,h,'LabelSpacing',100);
% % h.LineStyle = 'none';
% 
% cb = colorbar;
% ylabel(cb,'log_1_0(Viscosity, \eta [Pa s])');
% 
% set(gca, 'Layer','top'); box on; grid on;
% xlabel('Temperature, T [K]');
% ylabel('log_1_0(Grain Size, d [m])');
% 
% set(gca,'FontSize',fontsize);
% axis square
% 
% %%%%%%%%%%%%%%%%%
% % Shear Modulus Map
% %%%%%%%%%%%%%%%%%
% subplot(133)
% plot(T2(1,:),Gmod(1,:)/1e9,'linewidth',2);
% set(gca, 'Layer','top'); box on; grid on;
% xlabel('Temperature, T [K]');
% ylabel('Shear Modulus, \Gmod [GPa]');
% 
% set(gca,'FontSize',fontsize);
% axis square
% 
% 
% %%%%%%%%%%%%%%%%%
% % Grain size viscosity map
% %%%%%%%%%%%%%%%%%
% subplot(132)
% hold on;
% ylim([10,20]);
% [C,h] = contourf(T2,log10(eta),log10(d2),[-10:10]);
% clabel(C,h,'LabelSpacing',100);
% 
% cb = colorbar;
% % cb.TickLabels = labels;
% ylabel(cb,'log_1_0(Grain Size [m])');
% 
% set(gca, 'Layer','top'); box on; grid on;
% xlabel('Temperature, T [K]');
% ylabel('log_1_0(Viscosity, \eta [Pa s])');
% set(gca,'FontSize',fontsize);
% axis square
% 
% drawnow
% 
% for idx = 1 : length(h.TextPrims)
%     switch str2double(h.TextPrims(idx).String)
%         
%         case -7
%             h.TextPrims(idx).String=labels{1};
%         case -6
%             h.TextPrims(idx).String=labels{2};
%         case -5
%             h.TextPrims(idx).String=labels{3};
%         case -4
%             h.TextPrims(idx).String=labels{4};
%         case -3
%             h.TextPrims(idx).String=labels{5};
%         case -2
%             h.TextPrims(idx).String=labels{6};
%         case -1
%             h.TextPrims(idx).String=labels{7};
%         case -0
%             h.TextPrims(idx).String=labels{8};
%         case 1
%             h.TextPrims(idx).String=labels{9};
%             
%         otherwise
%     end
%     
% end
% 
% 
% %%%%%%%%%%%%%%%%%
% % Heat production controls
% %%%%%%%%%%%%%%%%%
% figure(2); clf; set(gcf,'color','w'); fontsize = 16;
% 
% %%%% IO %%%%
% subplot(131)
% contourf(T2,log10(d2),log10(q_e2_I),50);
% 
% cb = colorbar;
% % cb.TickLabels = labels;
% ylabel(cb,'log_1_0(Volumetric Dissipation Rate / \epsilon_0^2 [(W/m^3)/(m/m)^2])');
% 
% set(gca, 'Layer','top'); box on; grid on;
% ylabel('log_1_0(Grain Size [m])');
% xlabel('Temperature, T [K]');
% set(gca,'FontSize',fontsize);
% title('IO') 
% caxis([-5,5])
% axis square
% 
% %%%% EUROPA %%%%
% subplot(132)
% contourf(T2,log10(d2),log10(q_e2_E),50);
% 
% cb = colorbar;
% % cb.TickLabels = labels;
% ylabel(cb,'log_1_0(Volumetric Dissipation Rate / \epsilon_0^2 [(W/m^3)/(m/m)^2])');
% 
% set(gca, 'Layer','top'); box on; grid on;
% ylabel('log_1_0(Grain Size [m])');
% xlabel('Temperature, T [K]');
% set(gca,'FontSize',fontsize);
% title('EUROPA') 
% caxis([-5,5])
% axis square
% 
% 
% %%%% GANYMEDE %%%%
% subplot(133)
% contourf(T2,log10(d2),log10(q_e2_G),50);
% 
% cb = colorbar;
% % cb.TickLabels = labels;
% ylabel(cb,'log_1_0(Volumetric Dissipation Rate / \epsilon_0^2 [(W/m^3)/(m/m)^2])');
% 
% set(gca, 'Layer','top'); box on; grid on;
% ylabel('log_1_0(Grain Size [m])');
% xlabel('Temperature, T [K]');
% set(gca,'FontSize',fontsize);
% title('GANYMEDE') 
% caxis([-5,5])
% axis square
% 
% 
% 
% %%%%%%%%%%%%%%%%%
% % Melt Scaling
% %%%%%%%%%%%%%%%%%
% figure(3); clf; set(gcf,'color','w'); fontsize = 16;
% 
% h_m = linspace(0,1,1000); % Melted Depth
% V_m = linspace(0,1,1000); % Melted Volume
% R = 1; V0 = 1;
% qq0_h = ((R-h_m)/R).^(3);
% qq0_V = 1-V_m/V0;
% 
% hold on;
% plot(h_m,qq0_h);
% plot(V_m,qq0_V);
% 
% legend('h_m_e_l_t/R','V_m_e_l_t/V_0');
% 
% set(gca, 'Layer','top'); box on; grid on;
% set(gca,'FontSize',fontsize);
% xlabel('Scaling (See Label)')
% ylabel('$\dot{q}/\dot{q_0}$', 'Interpreter','latex','FontSize',fontsize+6)
% axis square
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
