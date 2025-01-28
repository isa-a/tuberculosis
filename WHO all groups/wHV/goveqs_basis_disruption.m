% Version to allow beta_factor inside this code, to simplify modelling
% of disruptions

function [out, lam] = goveqs_basis_disruption(t, in, M, notif_fac, beta_fac, i, prm, sel, agg)  

invec = in(1:i.nstates);

rHIV = interp1(0:length(prm.rHIV)-1, prm.rHIV, (t-1980));                  

% Normalise by populations
lam = M.lambda*invec/sum(invec)*beta_fac;
fac = interp1(2019+(0:length(notif_fac)), [1, notif_fac], t);
% Pull them together
allmat = M.lin + M.Dxlin*fac + rHIV*M.linHIV + lam*M.nlin;
out = allmat*invec;

% Implement deaths
morts = sum(M.mortvec,2).*invec;
out = out - morts;

% Implement births
births = sum(morts);
out(i.U.v0.h0) = out(i.U.v0.h0)+births;

% Get the auxiliaries
out(i.aux.inc)     = agg.inc*(sel.inc.*allmat)*invec;
out(i.aux.noti)    = agg.noti*(sel.noti.*allmat)*invec;
out(i.aux.mort(1)) = sum(M.mortvec(:,2).*invec);
out(i.aux.mort(2)) = sum(M.mortvec(:,3).*invec);
out(i.aux.Dx)      = agg.Dx*(sel.Dx.*allmat)*invec;