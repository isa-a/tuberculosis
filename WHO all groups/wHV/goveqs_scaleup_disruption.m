function out = goveqs_scaleup_disruption(t, in, M0, M1, times, notif_fac, beta_fac, i, prm, sel, agg)

scale = min(max((t-times(1))/(times(2)-times(1)),0),1);                    
Mt = M1; 
Mt.lin   = M0.lin + scale*(M1.lin-M0.lin);
Mt.Dxlin = M0.Dxlin + scale*(M1.Dxlin-M0.Dxlin);
out = goveqs_basis_disruption(t, in, Mt, notif_fac, beta_fac, i, prm, sel, agg);