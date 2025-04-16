clear all; load Model_setup; % load calibration_res_prev cov0;

obj  = @(x) get_objective3(x, ref, prm, gps, prm.contmat, rin_vec, lhd);

testvec = mean(prm.bounds);

[out, aux] = obj(testvec)
