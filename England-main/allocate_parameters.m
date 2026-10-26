function [p, r, prm] = allocate_parameters(x,p,r,xi,prm)

% x = x./scaling;

r.beta       = x(xi.r_beta);
r.migrout    = x(xi.r_migrout)/10;
p.migr_Lf    = x(xi.p_migr_Lf);
p.betadec    = x(xi.p_betadec)/100;
r.gamma_2015 = x(xi.r_gamma(1));
r.gamma_2020 = x(xi.r_gamma(1))*x(xi.r_gamma(2));        % NB: not currently used

r.progression  = r.progression0*x(xi.rf_progression)*[1, x(xi.rf_progr_migrVdom), 1];
r.reactivation = r.reactivation0*x(xi.rf_reactivation)*[1, x(xi.rf_progr_migrVdom), 1];
r.LTBI_stabil  = r.LTBI_stabil0*x(xi.rf_LTBIstabil);

r.muTB         = r.muTB0*x(xi.rf_muTB);
r.selfcure     = r.self_cure0*x(xi.rf_selfcure);

% p.migrmix = x(xi.migrmix);

prm.contmat = [1.0000, 0.5000; 0.5000, x(xi.migrmix)];

