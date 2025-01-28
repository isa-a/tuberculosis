function [r,p] = alloc_parameters(x,r,p,xi)

r.beta        = x(xi.r_beta(1))*[1 x(xi.r_beta(2))];
r.sym         = x(xi.r_sym);
r.cs          = x(xi.r_cs);
r.cs2         = x(xi.r_cs2);
p.Dx          = x(xi.p_Dx);
r.mort_TB     = r.mort_TB0*x(xi.rf_mort_TB);
r.self_cure   = r.self_cure0*x(xi.rf_self_cure);
% r.HIV_acqu = x(xi.r_HIV_acqu);
r.ART_init            = x(xi.r_ART_init);
r.HIV_mort            = x(xi.r_HIV_mort);

r.progression = r.progression0;
r.progression([2,3])  = r.progression0(1)*x(xi.p_HIV_relrate)*[1 0.4];

r.reactivation = r.reactivation0;
r.reactivation([2,3]) = r.reactivation0(1)*x(xi.p_HIV_relrate)*[1 0.4];

% if length(x) > xi.calib
%     r.self_cure([1,3]) = x(xi.r_self_cure);
% end