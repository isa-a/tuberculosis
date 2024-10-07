function [p,r] = allocate_parameters(x,p,r,xi)

r.beta           = x(xi.beta);
p.betadec        = x(xi.betadec);
r.gamma_2015     = x(xi.gamma(1))*[x(xi.p_relrate_gamma_chvad), 1];
r.gamma_2020     = x(xi.gamma(2))*[x(xi.p_relrate_gamma_chvad), 1];

tmp              = r.progression0*[1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
r.progression    = [tmp*x(xi.p_relrate_factor); tmp];

tmp              = r.reactivation0*[1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
r.reactivation   = [tmp*x(xi.p_relrate_factor); tmp];

p.LTBI_in_migrad = x(xi.p_LTBI_in_migrad);
p.LTBI_in_migrch = x(xi.p_LTBI_in_migrad)/x(xi.p_relLTBI_inmigr_advch);

% r.vuln           = x(xi.r_vuln);
% p.relbeta_vuln   = 1;

r.ageing         = x(xi.ageing);

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
%r.vuln           = x(xi.r_vuln);
%p.relbeta_vuln   = 1;