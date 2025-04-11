function [p,r] = allocate_parameters(x,p,r,xi,scaling)

% x = x./scaling;

r.beta           = x(xi.beta);

r.migrout        = x(xi.r_migrout);

% p.betadec        = x(xi.betadec);
p.betadec        = 0;

r.gamma_2015     = x(xi.gamma(1))*[x(xi.p_relrate_gamma_chvad), 1];
r.gamma_2020     = x(xi.gamma(2))*[x(xi.p_relrate_gamma_chvad), 1];        % NB: not currently used

% tmp              = r.progression0*[1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
tmp              = r.progression0*[1, x(xi.p_relrate(1)), 1];
r.progression    = [tmp*x(xi.p_relrate_factor); tmp];

% tmp              = r.reactivation0*[1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
tmp              = r.reactivation0*[1, x(xi.p_relrate(1)), 1];
r.reactivation   = [tmp*x(xi.p_relrate_factor); tmp];

p.LTBI_in_migrad = x(xi.p_LTBI_in_migrad);
p.LTBI_in_migrch = x(xi.p_LTBI_in_migrad)/x(xi.p_relLTBI_inmigr_advch);

% r.vuln           = x(xi.r_vuln_sc)/500;
% p.relbeta_vuln   = 1;                        % NB: not currently used in make_model

r.ageing         = x(xi.r_ageing_sc)/10;

r.muTB           = r.muTB*x(xi.mort_factor);
%r.ch_mort        = x(xi.ch_mort);

% r.progression = [
%     r.progression0 * [1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))]; % ch
%     r.progression0 * p_relrate_factor * [1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))]  % ad
% ];
% 
% r.reactivation = [
%     r.reactivation0 * [1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))]; % ch
%     r.reactivation0 * p_relrate_factor * [1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))]  % ad
% ];

%r.migr           = x(xi.r_migr);