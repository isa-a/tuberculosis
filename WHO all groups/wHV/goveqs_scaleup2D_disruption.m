function out = goveqs_scaleup2D_disruption(t, in, M0, M1, M2, times, notif_fac, i, prm, sel, agg)

scale = max((t-times(:,1))./(times(:,2)-times(:,1)),0);                    
scale(2) = min(scale(2),1);                                                % Only applies on second argument, so that model can continue extrapolating ART coverage into future
Mt = M1; 
Mt.lin   = M0.lin + scale(1)*(M1.lin-M0.lin) + scale(2)*(M2.lin-M0.lin);
Mt.Dxlin = M0.Dxlin + scale(1)*(M1.Dxlin-M0.Dxlin) + scale(2)*(M2.Dxlin-M0.Dxlin);
% out = goveqs_basis2(t, in, Mt, i, s, p, sel, agg);
out = goveqs_basis_disruption(t, in, Mt, notif_fac, i, prm, sel, agg);