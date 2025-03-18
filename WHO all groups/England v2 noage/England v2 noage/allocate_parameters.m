function [p,r] = allocate_parameters(x,p,r,xi,scaling)

% x = x./scaling;

r.beta           = x(xi.beta);
p.betadec        = x(xi.betadec);
% r.gamma_2015     = x(xi.gamma(1))*[x(xi.p_relrate_gamma_chvad), 1];
% r.gamma_2020     = x(xi.gamma(2))*[x(xi.p_relrate_gamma_chvad), 1];

r.gamma_2015     = x(xi.gamma(1));
r.gamma_2020     = x(xi.gamma(2));

tmp              = r.progression0*[1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
r.progression    = tmp;
% r.progression    = [tmp*x(xi.p_relrate_factor); tmp];

tmp              = r.reactivation0*[1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
r.reactivation   = tmp;
% r.reactivation   = [tmp*x(xi.p_relrate_factor); tmp];

p.LTBI_in_migrad = x(xi.p_LTBI_in_migrad);
% p.LTBI_in_migrch = x(xi.p_LTBI_in_migrad)/x(xi.p_relLTBI_inmigr_advch);

r.vuln           = x(xi.r_vuln_sc)/500;
p.relbeta_vuln   = 1;

% r.ageing         = x(xi.r_ageing_sc)/10;