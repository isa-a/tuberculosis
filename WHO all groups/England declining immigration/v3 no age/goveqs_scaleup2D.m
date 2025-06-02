function out = goveqs_scaleup2D(t, in, M0, M1, M2, rin_vec, times, i, s, p, r, prm, sel, agg, equil)

    % ─────────────────────────────────────────────────────────────────
    % 1) Override LTBI_in_migrad to 0.35 once t ≥ 2020:
    if t >= 2020
        p.LTBI_in_migrad = 0.6;
    end
    % ─────────────────────────────────────────────────────────────────

    scale = max((t - times(:,1)) ./ (times(:,2) - times(:,1)), 0);  
    % <--- NEW: removed min(,1), so that model can continue extrapolating ART coverage into future
    scale(1) = min(scale(1), 1);

    Mt = M1;
    Mt.lin = M0.lin + scale(1)*(M1.lin - M0.lin) + scale(2)*(M2.lin - M0.lin);

    out = goveqs_basis3(t, in, i, s, Mt, rin_vec, agg, sel, r, p, equil);
end
