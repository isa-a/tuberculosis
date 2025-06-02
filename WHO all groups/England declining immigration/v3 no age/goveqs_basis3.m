function out = goveqs_basis2(t, in, i, s, M, rin_vec, agg, sel, r, p, equil)

    out = zeros(length(in),1);
    invec = in(1:i.nstates);

    % Prepare population denominators
    tmp = M.denvec * invec;
    den = sum(tmp .* M.denvec, 1)';
    den(den==0) = Inf;

    % New infections
    lam = M.lam * (invec ./ den);

    % ─────────────────────────────────────────────────────────────────
    % ONLY CHANGE: use 0.35 for LTBI_in_migrad if t ≥ 2020, else use calibrated p.LTBI_in_migrad
    allmat = M.lin ...
           + M.LTBIlin * (((t >= 2020) * 0.35 + (t < 2020) * p.LTBI_in_migrad)) ...
               * (1 - p.LTBIdec)^(max((t - 2010), 0)) ...
           + lam(1) * M.nlin.ad.dom.ds       ...
           + lam(2) * M.nlin.ad.migr_rect.ds ...
           + lam(2) * M.nlin.ad.migr_long.ds;
    % ─────────────────────────────────────────────────────────────────

    out(1:i.nstates) = allmat * invec;

    % Mortality
    morts = M.mort .* invec;
    out(1:i.nstates) = out(1:i.nstates) - sum(morts, 2);

    % Births into UK population
    dom_morts = sum(sum(morts([s.dom], :)));
    out(i.U.ad.dom) = out(i.U.ad.dom) + dom_morts;

    if equil
        r_in = rin_vec(1);
    else
        r_in = rin_vec(min(max(ceil(t - 2019), 1), length(rin_vec)));
    end

    % Migration into UK (unchanged)
    out(1:i.nstates) = out(1:i.nstates) + r_in * M.migrentries;

    if sum(invec) > 5
        keyboard;
    end

    % Auxiliaries
    out(i.aux.inc)        = agg.inc * (sel.inc .* allmat) * invec;
    out(i.aux.mort)       = sum(morts(:, 2));
    out(i.aux.nTPT)       = sum((sel.nTPT .* allmat) * invec);
    out(i.aux.ch_notifs)  = sum((sel.ch_notifs .* allmat) * invec);
end
