function [p,r,prm] = allocate_parameters(x,p,r,prm,xi)

r.beta           = x(xi.beta);
p.betadec        = x(xi.betadec);
r.gamma_2015     = x(xi.gamma(1))*[x(xi.p_relrate_gamma_chvad), 1];
r.gamma_2020     = x(xi.gamma(2))*[x(xi.p_relrate_gamma_chvad), 1];

r.progression = r.progression0 * ones(2,2,3); 
r.progression(:,:,2:3) = r.progression0 * x(xi.p_HIV_relrate) .* repmat(reshape([1 0.4], 1, 1, []), 2, 2, 1);

r.reactivation = r.reactivation0 * ones(2,2,3); 
r.reactivation(:,:,2:3) = r.reactivation0 * x(xi.p_HIV_relrate) .* repmat(reshape([1 0.4], 1, 1, []), 2, 2, 1);


% r.vuln           = x(xi.r_vuln_sc)/500;
% p.relbeta_vuln   = 1;
r.ageing         = x(xi.r_ageing_sc);

tmp_c = prm.contmat_age;
prm.contmat_age = [tmp_c(1, :) * x(xi.contmat_factor); tmp_c(2, :)];
prm.contmat = zeros(4, 4); % Adjust size if needed
for age_row = 1:2 % rows in age
    for age_col = 1:2 % cols in age
        for born_row = 1:3 % rows in born
            for born_col = 1:3 % cols in born
                % Calculate position in combined matrix
                row = (born_row - 1) * 2 + age_row;
                col = (born_col - 1) * 2 + age_col;
                % Multiply
                prm.contmat(row, col) = prm.contmat_born(born_row, born_col) * prm.contmat_age(age_row, age_col);
            end
        end
    end
end


r.ART_init            = x(xi.r_ART_init);
r.HIV_mort            = x(xi.r_HIV_mort);