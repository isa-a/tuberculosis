function [out, lam] = goveqs_basis(t, in, M, i, s, p, sel, agg)

invec = in(1:i.nstates);

% Normalise by populations
lam = M.lambda*invec/sum(invec);
allmat = M.lin + lam(1)*M.nlin.DS + lam(2)*M.nlin.MDR;

out = allmat*invec;

% Implement deaths
morts = M.mortvec.*invec;
out = out - morts;

% Implement births
births = sum(morts);
out(i.U) = out(i.U)+births;

% Get the auxiliaries
out(i.aux.inc(1:2)) = agg.inc*(sel.inc.*allmat)*invec;
out(i.aux.inc(3))   = sum((sel.acqu.*allmat)*invec);
