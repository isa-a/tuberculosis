function [p,r,prm] = allocate_parameters(x,p,r,prm,xi)

r.beta           = x(xi.beta);
p.betadec        = x(xi.betadec);
r.gamma_2015     = x(xi.gamma(1))*[x(xi.p_relrate_gamma_chvad), 1];
r.gamma_2020     = x(xi.gamma(2))*[x(xi.p_relrate_gamma_chvad), 1];

% tmp              = r.progression0*[1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
% r.progression    = [tmp*x(xi.p_relrate_factor); tmp];

tmp = r.progression0 * [1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
tmp2 = [tmp * x(xi.p_relrate_factor); tmp];
hiv_factors = [1, x(xi.HIVfactor), x(xi.HIVfactor) * 0.4];
hiv_factors = reshape(hiv_factors, [1, 1, numel(hiv_factors)]);
r.progression = tmp2 .* hiv_factors;

% tmp              = r.reactivation0*[1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
%r.reactivation   = [tmp*x(xi.p_relrate_factor); tmp];

tmp = r.reactivation0 * [1, x(xi.p_relrate(1)), 1, x(xi.p_relrate(2))];
tmp2 = [tmp * x(xi.p_relrate_factor); tmp];
hiv_factors = [1, x(xi.HIVfactor), x(xi.HIVfactor) * 0.4];
hiv_factors = reshape(hiv_factors, [1, 1, numel(hiv_factors)]);
r.reactivation = tmp2 .* hiv_factors;

% p.LTBI_in_migrad = x(xi.p_LTBI_in_migrad);
% p.LTBI_in_migrch = x(xi.p_LTBI_in_migrad)/x(xi.p_relLTBI_inmigr_advch);

r.vuln           = x(xi.r_vuln_sc)/500;
p.relbeta_vuln   = 1;

r.ageing         = x(xi.r_ageing_sc)/10;
%r.ch_mort        = x(xi.ch_mort);

tmp_c = prm.contmat_age;
prm.contmat_age = [tmp_c(1, :) * x(xi.contmat_factor); tmp_c(2, :)];

prm.contmat      = zeros(6, 6);
% go through each element
for age_row = 1:2                                                           % rows in age
    for age_col = 1:2                                                       % cols in age
        for born_row = 1:2                                                  % rows in born
            for born_col = 1:2                                              % cols in born
                % calc position in combined matrix
                row = (born_row-1)*2 + age_row;                             % correctly scale current rows into new matrix
                col = (born_col-1)*2 + age_col;                             % correctly scale current cols into new matrix
                % multiply
                prm.contmat(row, col) = prm.contmat_born(born_row, born_col) * prm.contmat_age(age_row, age_col);
            end
        end
    end
end

r.HIVincdpeak         = x(xi.HIVincdpeak);
r.HIVincdnow         = x(xi.HIVincdnow);

r.ARTnow         = x(xi.r_ARTnow);

r.muHIV          = x(xi.muHIV);
r.muTBHIV        = x(xi.muTBHIV);

r.HIV        = x(xi.HIV);
r.ART        = x(xi.ART);

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
