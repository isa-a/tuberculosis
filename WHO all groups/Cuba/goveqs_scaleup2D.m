function out = goveqs_scaleup2D(t, in, M1, M2, M3, M4, M5, times, i, s, p, r, prm, sel, agg)

scale = max((t-times(:,1))./(times(:,2)-times(:,1)),0);                    % <--- NEW: removed min(,1), so that model can continue extrapolating ART coverage into future
scale(1) = min(scale(1),1);

Mt = M2; Mt.lin = M2.lin + scale(1)*(M2.lin-M1.lin) + scale(2)*(M2.lin-M1.lin) + scale(3)*(M3.lin-M1.lin) + scale(4)*(M4.lin-M1.lin) + scale(5)*(M5.lin-M1.lin);
out = goveqs_basis3(t, in, i, s, Mt, agg, sel, r, p);