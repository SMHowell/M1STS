%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M,IN] = initializeGrid(IN,BOD)
%%%%%%%%%%%%%%%%%%%%%%
% Interfaces
%%%%%%%%%%%%%%%%%%%%%%
M.rOcn = BOD.R;  % Outer radius of ocean surface
M.rSil = BOD.R;  % Outer radius of seafloor
M.rIrn = IN.H0_irn; % Outer radius of core


%%%%%%%%%%%%%%%%%%%%%%
% Differentiation-related Geometry
%%%%%%%%%%%%%%%%%%%%%%
% Find largest plausible H2O layer to ensure well-resolved
if isempty(IN.Hmax_H2O)
    rSil_min    = floor((BOD.R^3 - (3/(4*pi))*(IN.fm0_H2O*BOD.m)/BOD.rhoIce_0)^(1/3));
    IN.Hmax_H2O = BOD.R - rSil_min;
else
    rSil_min   = BOD.R - IN.Hmax_H2O;
    IN.fm0_H2O = (4/3)*pi*(BOD.R^3-rSil_min^3)*BOD.rhoIce_0/BOD.m; 
end

% Find largest plausible iron layer to ensure well-resolved
if isempty(IN.Hmax_irn)
    rIrn_max    = ((3/(4*pi))*(IN.fm0_irn*BOD.m)/BOD.rhoIrn_0)^(1/3);
    IN.Hmax_irn = rIrn_max;
else
    rIrn_max    = IN.Hmax_irn;
    IN.fm0_irn  = (4/3)*pi*IN.Hmax_irn^3*BOD.rhoIrn_0/BOD.m; 
end


%%%%%%%%%%%%%%%%%%%%%%
% Model Grid
%%%%%%%%%%%%%%%%%%%%%%
r1   = rSil_min; r2 = BOD.R; dz = min(IN.dz_H2O,IN.dz_sil);
rH2O = linspace(r1, r2, ceil((r2-r1)/dz));

r1   = rIrn_max; r2 = rSil_min; dz = IN.dz_sil;
rSil = linspace(r1, r2, ceil((r2-r1)/dz));

r1   = 0; r2 = rIrn_max; dz = min(IN.dz_irn,IN.dz_sil);
rIrn = linspace(r1, r2, ceil((r2-r1)/dz));

% Radial array with space ghost node 1 km off the surface
M.r   = sort(unique([rH2O,rSil,rIrn,BOD.R+IN.dz_H2O]));
M.dr  = diff(M.r);    % Element height [m]
M.z   = BOD.R - M.r;  % Depth [m]

M.Nz = numel(M.r);    % Number of nodes in vertical direction
M.ind = 1:M.Nz-1;     % A useful index for referencing values 

% Define the location of elements in thermal calculations. Note that M.r_s
% converges to the midpoint between nodes when ice thickness << body radius
M.r_s    = sqrt(M.r(M.ind) .* M.r(M.ind+1));  % Radius from center of body [m]
M.r_s(1) = sqrt(M.dr(1));                     % Arbitrary correction     
M.z_s    = BOD.R - M.r_s;                     % Depth [m]

% Define element volumes, surface areas, and length scales
M.V     = zeros(1,M.Nz); % Representative volumes on nodes
M.V(2:end-1) = (4/3) * pi * (M.r_s(2:end).^3-M.r_s(1:end-1).^3);
M.V(1)  = (4/3) * pi * M.r_s(1).^3;
M.V(end)= (4/3) * pi * (M.r(end).^3-M.r_s(end).^3);

M.V_s = (4/3) * pi * (M.r(M.ind+1).^3-M.r(M.ind).^3);
M.A   = 4 * pi * M.r.^2;
M.L   = M.V_s./M.A(M.ind+1);  % characteristic length scale

% Track element containing interface
M.iOcnTop = find((M.rOcn - M.r)>=0,1,'last'); % Ocean top interface element index
M.iOcnBot = find((M.rSil - M.r)>=0,1,'last'); % Ocean bottom interface element index
M.iSilBot = find((M.rIrn - M.r)>=0,1,'last'); % Ocean bottom interface element index

%%%%%%%%%%%%%%%%%%%%%%
% The volumetric weightings are pre-computed to save time in later
% interpolations
%%%%%%%%%%%%%%%%%%%%%%
ind     = M.ind; % Indices for elements b/w nodes
M.fV1   = M.V(ind)  ./(M.V(ind+1)+M.V(ind)); % Weighting for node i
M.fV2   = M.V(ind+1)./(M.V(ind+1)+M.V(ind)); % Weighting for node i+1

ind     = 2:M.Nz-1; % Indices for nodes b/w elements
M.fV1_s = M.V_s(ind-1)./(M.V_s(ind-1) + M.V_s(ind)); % Weighting for element i+1
M.fV2_s = M.V_s(ind)  ./(M.V_s(ind-1) + M.V_s(ind)); % Weighting for element i+1

end














