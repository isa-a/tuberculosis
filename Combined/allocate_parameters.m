function [p,r] = allocate_parameters(x,p,r,xi)

r.beta         = x(xi.beta);
r.gamma        = x(xi.gamma);
p.birth        = x(xi.p_birth);
p.kLf          = x(xi.p_kLf);
r.adu_mort     = x(xi.adu_mort);
r.progression  = r.progression0*[1, x(xi.eld_relrisk)];
r.reactivation = r.reactivation0*[1, x(xi.eld_relrisk)];