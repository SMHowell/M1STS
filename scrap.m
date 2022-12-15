syms g1 g2 g3 f1 f2 f3 rho1 rho2 rho3 

eqn1 = g1 == f1 / ((1-f1)*(g2*rho2 + g3*rho3)/rho1 + f1);
eqn2 = g2 == f2 / ((1-f2)*(g1*rho1 + g3*rho3)/rho2 + f2);
eqn3 = g3 == f3 / ((1-f3)*(g1*rho1 + g2*rho2)/rho3 + f3);

sol = solve([eqn1, eqn2, eqn3], [g1, g2, g3])





% Apply to volumes
% Get ice from rock-iron
fm_RI_irn = fm_irn./(fm_irn + fm_sil);
fV_RI_irn = fm_RI_irn./(fm_RI_irn+(1-fm_RI_irn).*rho_irn./rhoSil_s);
rho_RI    = fV_RI_irn.*rho_irn + (1-fV_RI_irn).*rhoSil_s;

% Now calculate ice volume fraction!
fV_H2O = fm_H2O./(fm_H2O+(1-fm_H2O).*rhoH2O_s./rho_RI);
fV_H2O(isnan(fV_H2O)) = 1;

% Rock-Ice portion
fm_RH_H2O = fm_H2O./(fm_H2O + fm_sil);
fV_RH_H2O = fm_RH_H2O./(fm_RH_H2O+(1-fm_RH_H2O).*rhoH2O_s./rhoSil_s);
rho_RH    = fV_RH_H2O.*rhoH2O_s + (1-fV_RH_H2O).*rhoSil_s;

% Now calculate iron volume fraction!
fV_irn = fm_irn./(fm_irn+(1-fm_irn).*rho_irn./rho_RH);
fV_irn(isnan(fV_irn)) = 1;

% Ice-Iron portion
fm_IH_H2O = fm_H2O./(fm_H2O + fm_irn);
fV_IH_H2O = fm_IH_H2O./(fm_IH_H2O+(1-fm_IH_H2O).*rhoH2O_s./rho_irn);
rho_IH    = fV_IH_H2O.*rhoH2O_s + (1-fV_IH_H2O).*rho_irn;

% Now calculate silicate volume fraction!
fV_sil = fm_sil./(fm_sil+(1-fm_sil).*rhoSil_s./rho_IH);
fV_sil(isnan(fV_sil)) = 1;