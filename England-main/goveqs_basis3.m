function out = goveqs_basis2(t, in, i, s, M, rin_vec, agg, sel, r, p, equil)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% New infections
lam = M.lam*invec*(1-p.betadec*(1-equil))^(max((t-2010),0));
% lam = M.lam*(invec./sum(invec))*(1-p.betadec)^(max((t-2010),0));

allmat = M.lin + M.LTBIlin + ...
        lam(1)*M.nlin.dom + ...
        lam(2)*M.nlin.migr_rect + ...
        lam(2)*M.nlin.migr_long;

out(1:i.nstates) = allmat*invec;

% Mortality
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);

% Replenishing UK population
dom_morts = sum(sum(morts(s.dom,:)));
out(i.U.dom) = out(i.U.dom) + dom_morts;

% Replenishing migrant population
if equil
    r_in = rin_vec(1);
else
    r_in = rin_vec(min(max(round(t) - 2019 + 1,1),length(rin_vec)));
end
out(1:i.nstates) = out(1:i.nstates) + r_in*M.migrentries;

% if sum(invec)>5
%     keyboard;
% end


% Auxiliaries
out(i.aux.inc)        = agg.inc*(sel.inc.*allmat)*invec;
tmp1                  = agg.incsources*((sel.L2I.*allmat)*invec);
tmp2                  = agg.incsources*((sel.P2I.*allmat)*invec);
tmp3                  = agg.incsources*((sel.R2I.*allmat)*invec);
tmp4                  = agg.incsources*((sel.T2I.*allmat)*invec);
out(i.aux.incsources) = [tmp1; tmp2; tmp3; tmp4];
out(i.aux.mort)       = sum(morts(:,2));
out(i.aux.nTPT)       = sum((sel.nTPT.*allmat)*invec);
