%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiphase 1-dimensional Satellite Thermal Solver
% Spherically symmetric solver for the heat diffusion equation, explicit
% finite difference in time. Tailored for applications to planetary
% evolution with considerations for phase change, reservoir freezing, and 
% eruption.
% Sam Howell, Elodie Lesage, Julia Miller
% (C)2022 California Institute of Technology. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = getPressure(M,BOD,MAT)

% Guestimate -- *** REPLACE THIS ***
P_1 = (BOD.R-M.rSil)  * MAT.H2O.m.rho0 * BOD.g;
P_2 = P_1 + (M.rSil-M.rIrn) * MAT.SIL.s.rho0 * BOD.g;
P_3 = P_2 + M.rIrn * MAT.IRN.s.rho0 * BOD.g;

indH2O = find(M.r<=BOD.R  & M.r>M.rSil);
indSil = find(M.r<=M.rSil & M.r>M.rIrn);
indIrn = find(M.r<=M.rIrn);

P = zeros(1,M.Nz);
P(indH2O) = linspace(P_1,  0,numel(indH2O)); 
P(indSil) = linspace(P_2,P_1,numel(indSil)); 
P(indIrn) = linspace(P_3,P_2,numel(indIrn)); 

M.P = P;

end









































