clear all; load optim_res1000.mat;

% rin_vec = 0.01*ones(size(rin_vec));

% rin_vec(end-2:end) = rin_vec(end-3);
obj = @(x) get_objective3(x, ref, prm, gps, prm.contmat, rin_vec, lhd);

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

inds = find(outs==max(outs));
xx   = xsto(inds(1),:);

% --- Get the past incidence
[out,aux] = obj(xx);

init  = aux.init;
incd  = aux.incd;
incdt = aux.incdt;

% --- Get the future incidence
[p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);

% r0.gamma = r0.gamma_2020;
p0.prev_in_migr = 0;
r0.gamma        = r0.gamma_2015;
r0.TPT          = [0 r0.TPT2020rec 0];
M0              = make_model(p0, r0, i, s, gps, prm0.contmat);

geq = @(t,in) goveqs_basis3(t, in, i, s, M0, rin_vec, agg, sel, r0, p0, false);

[t, soln] = ode15s(geq, 2021:2041, init, opts);
sdiff  = diff(soln, [], 1);
sfin   = soln(end,:);
pops   = sum(soln(:,1:i.nstates),2);
incsto = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);

incd
incsto(1)

return;

midpt = false; 
if midpt
    xs = x0sto(2,:);
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
      
    init    = aux.soln(end, :);
    incd(ii) = aux.incd;

    [p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);
    % r0.gamma = r0.gamma_2020;
    % p0.prev_in_migr = 0;
    r0.gamma        = r0.gamma_2015;
    % r0.TPT          = [0 r0.TPT2020rec 0];
    M0              = make_model(p0,r0,i,s,gps,prm0.contmat);

    geq = @(t,in) goveqs_basis3(t, in, i, s, M0, rin_vec, agg, sel, r0, p0, false);

    [t, soln] = ode15s(geq, 2022:2041, init, opts);
    sdiff = diff(soln, [], 1);
    sfin  = soln(end,:);
    pops  = sum(soln(:,1:i.nstates),2);
    incsto(:, ii) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
    
end
fprintf('\n');

mat = prctile(incsto,[2.5,50,97.5],2)

prctile(incd,[2.5,50,97.5])

figure; plot(mat)

