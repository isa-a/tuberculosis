function out = goveqs_scaleup1D(t, in, M0, M1, times, i, s, p, r, prm, sel, agg)

scale = max((t-times(:,1))./(times(:,2)-times(:,1)),0);                    % <--- NEW: removed min(,1), so that model can continue extrapolating ART coverage into future
scale(1) = min(scale(1),1);

Mt = M1; Mt.lin = M0.lin + scale(1)*(M1.lin-M0.lin);
out = goveqs_basis3(t, in, i, s, Mt, agg, sel, r, p);