clear all; load tmp26.mat; load Model_setup.mat;
% rin_vec  = [rin_vec, rin_vec(end)*5/9];
saving = 1;
obj = @(x) get_objective3(x, ref, prm, gps, prm.contmat, rin_vec, lhd);
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
midpt = false;
if midpt
    xs = x3;
else
    ix0 = size(xsto,1)/2;
    nx  = 20;
    dx  = round(ix0/nx);
    xs  = xsto(ix0:dx:end,:);
end
mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    xx = xs(ii,:);
    [out,aux] = obj(xx);
    init = aux.init;
    init = aux.soln(2010:2022==2014,:);
    [p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);
    p0.prev_in_migr = 0;
    r0.gamma = r0.gamma_2015;
    r0.TPT = [0 r0.TPT2020rec 0];
    M0 = make_model(p0, r0, i, s, gps, prm0.contmat);
% Improved Tx outcomes
    r1 = r0; p1 = p0;
% Find new treatment outcomes
    vec = [r0.Tx, r0.ltfu, r0.muTx]; props = vec/sum(vec);
    props(end) = 0.01; props(2) = props(2)/2; props(1) = 1 - sum(props(2:3));
    newrates = r0.Tx*props/props(1);
    r1.ltfu = newrates(2);
    r1.muTx = newrates(3);
    M1 = make_model(p1, r1, i, s, gps, prm0.contmat);
% Enhanced TPT, recent migrants
    r2 = r1; p2 = p1;
    %r2.TPT = -log(1-0.25) * [0 1 0 0];
    r2.TPT = -log(1-0.25) * [0 1 0];
    M2 = make_model(p2, r2, i, s, gps, prm0.contmat);
% ACF in foreign-born
    r3 = r2; p3 = p2;
    r3.ACF = -log(1-0.25) * [0 1 1 1];
    M3 = make_model(p3, r3, i, s, gps, prm0.contmat);
% ACF in domestic and foreign-born
    r4 = r3; p4 = p3;
    r4.ACF = -log(1-0.25) * [1 1 1 1];
    M4 = make_model(p4, r4, i, s, gps, prm0.contmat);
% Pre-entry TPT
    r5 = r4; p5 = p4;
    p5.migrTPT = 0.6;
    M5 = make_model(p5, r5, i, s, gps, prm0.contmat);
    models = {M0, M1, M2, M3, M4, M5};
% prev_soln = init;
    M0_end_soln = [];
for mi = 1:length(models)
        geq = @(t,in) goveqs_scaleupb(t, in, i, s, M0, models{mi}, rin_vec, [2027 2035], agg, prm0, sel, r0, p0, false);
        [t, soln] = ode15s(geq, 2014:2041, init, opts);
        sdiff = diff(soln, [], 1);
        pops = sum(soln(:,1:i.nstates),2);
        incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
        mrtsto(:, ii, mi) = sdiff(:, i.aux.mort) * 1e5./pops(1:end-1);
        prvsto(:, ii, mi) = sum(soln(:,s.allI),2) * 1e5./pops;
if mi == 6
            incsource(ii,:) = sdiff(end,i.aux.incsources);
end
end
end
fprintf('\n');

if saving 
    save intvn_res2026_04;
end