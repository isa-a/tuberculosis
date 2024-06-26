function out = goveqs_basis2(t, in, i, s, M, agg, sel, r, p)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% Prepare population denominators
tmp = M.denvec*invec;
den = sum(tmp.*M.denvec,1)';
den(den==0) = Inf;

% New infections
lam = M.lam*(invec./den)*(1-p.betadec)^(max((t-2010),0));
% lam = M.lam*(invec./sum(invec))*(1-p.betadec)^(max((t-2010),0));

% Full model
allmat = M.lin + ...
        lam(1)*M.nlin.ch.dom.ds             + lam(2)*M.nlin.ch.dom.rr + ...
        lam(3)*M.nlin.ad.dom.ds             + lam(4)*M.nlin.ad.dom.rr + ...
        lam(5)*M.nlin.ch.migr_rect.ds       + lam(6)*M.nlin.ch.migr_rect.rr + ...
        lam(7)*M.nlin.ch.migr_long.ds       + lam(8)*M.nlin.ch.migr_long.rr + ...
        lam(9)*M.nlin.ad.migr_rect.ds       + lam(10)*M.nlin.ad.migr_rect.rr + ...
        lam(11)*M.nlin.ad.migr_long.ds      + lam(12)*M.nlin.ad.migr_long.rr + ...
        lam(13)*M.nlin.ch.vuln.ds           + lam(14)*M.nlin.ch.vuln.rr + ...
        lam(15)*M.nlin.ad.vuln.ds           + lam(16)*M.nlin.ad.vuln.rr;
out(1:i.nstates) = allmat*invec;

% Mortality
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);

% Births into UK population
dom_morts = sum(sum(morts([s.dom,s.vuln],:)));
out(i.U.dom) = out(i.U.dom) + dom_morts;

% Migration out of UK
out(s.migr) = out(s.migr) - r.migr*invec(s.migr)/sum(invec(s.migr));

% Migration into UK
inmigr = sum(sum(morts(s.migr,:))) + r.migr;
% vec = [1-p.LTBI_in_migr, (1-p.migrTPT)*p.LTBI_in_migr*[0.02 0.98], p.migrTPT*p.LTBI_in_migr*[0.02 0.98]]';
% out(s.migrstates) = out(s.migrstates) + inmigr*vec;
out(1:i.nstates) = out(1:i.nstates) + inmigr.*M.migrentries;

% % Migration
% out(1:i.nstates) = out(1:i.nstates) + M.migration;

% Auxiliaries
out(i.aux.inc)        = agg.inc*(sel.inc.*allmat)*invec;
tmp1                  = agg.incsources*((sel.L2I.*allmat)*invec);
tmp2                  = agg.incsources*((sel.P2I.*allmat)*invec);
tmp3                  = agg.incsources*((sel.R2I.*allmat)*invec);
tmp4                  = agg.incsources*((sel.T2I.*allmat)*invec);
out(i.aux.incsources) = [tmp1; tmp2; tmp3; tmp4];
out(i.aux.mort)       = sum(morts(:,2));
out(i.aux.nTPT)       = sum((sel.nTPT.*allmat)*invec);