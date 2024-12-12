clear all; load Model_setup; % load calibration_res_prev cov0;

obj  = @(x) get_objective2(x, ref, prm, gps, prm.contmat, lhd);
nobj = @(x) -obj(x);

nsam = 100; 
xsam = repmat(prm.bounds(1,:),nsam,1) + diff(prm.bounds).*lhsdesign(nsam,size(prm.bounds,2));

% obj(xsam(1,:));

% xx = xsam(1,:);
% [p,r] = allocate_parameters(xx, p, r, xi);
% init = zeros(1,i.nx);
% seed = 1e-5;
% init(i.U.dom)       = 1 - 0.168 - seed;
% init(i.U.migr_rect) = 0.168;
% init(i.I.dom)       = Setupseed;
% 
% p0 = p; r0 = r; p0.betadec = 0;
% M0 = make_model(p0, r0, i, s, gps);
% geq0 = @(t,in) goveqs_basis3(t, in, i, s, M0, agg, sel, r0, p0);
% [t0, soln0] = ode15s(geq0, [0:1e4], init, odeset('NonNegative',1:5));

mk = round(nsam/25);
for ii = 1:nsam
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    [outs(ii),~,msg(ii)] = obj(xsam(ii,:));
end
% Order by fit
mat  = sortrows([outs; 1:nsam]',-1);
ord  = mat(:,2);
xord = xsam(ord,:);

options = optimset(PlotFcn=@optimplotfval);
for ii = 1:10
    tmp = fminsearch(nobj,xord(ii,:),options);
    x0sto(ii,:) = fminsearch(nobj,tmp,options);
end

% x1 = fminsearch(nobj,x0,options);
% x2 = fminsearch(nobj,x1,options);
% save optim_res_noVULN5_v2;

save optim_res_MAIN
% [xsto, outsto] = MCMC_adaptive2(obj, x0sto(2,:), 1000, 1, [], true);


return;


xopt = zeros(3, size(xord, 2));
fvals_opt = zeros(1, 3);    


for ii = 1:3
    options = optimset('PlotFcn',@optimplotfval);
    [xopt(ii,:), fval] = fminsearch(nobj, xord(ii,:), options);
    fvals_opt(ii) = -fval;  % Store the maximized objective function value
end


% return;

options = optimset(PlotFcn=@optimplotfval);
x0 = fminsearch(nobj,xord(1,:),options);
x1 = fminsearch(nobj,x0,options);
x2 = fminsearch(nobj,x1,options);
x3 = fminsearch(nobj,x2,options);



% Perform MCMC
[xsto, outsto] = MCMC_adaptive33(obj, x0, 1e2, 1, [], [], [], 1);

inds = find(outsto==max(outsto));
x_new = xsto(inds(1),:);


cov0 = cov(xsto);
[xsto, outsto] = MCMC_adaptive33(obj, x_new, 1e4, 1, [], [], cov0, 1);


nreps = 4;
niter = [1, 1, 1, 5]*2e3;
for ii = 1:nreps
    [xsto, outsto] = MCMC_adaptive2(obj, x2, niter(ii), 1, [], 1);
    inds = find(outsto==max(outsto));
    x2 = xsto(inds(1),:);
    cov0 = cov(xsto);
    fprintf('\n');
end


save calibration_res_december;
return;

x0 = fminsearch(nobj,xord(1,:),optimset('Display','iter'));
x0 = fminsearch(nobj,x0,optimset('Display','iter'));
x0 = fminsearch(nobj,x0,optimset('Display','iter'));

x0_init = xord;
save calibration_res_isa_new;

return;

% % Perform MCMC
% [xsto, outsto] = MCMC_adaptive(obj, x0, 5e4, 1, [], [], cov0, 1);
% 
% inds = find(outsto==max(outsto));
% x0 = xsto(inds(1),:);
% 
% [out, aux] = obj(x0);
% sfin = aux.soln(end,:);
% sum(sfin(intersect(s.migr,[s.Lf,s.Ls])))/sum(sfin(s.migr))
% 
% 
% cov0 = cov(xsto);
% [xsto, outsto] = MCMC_adaptive(obj, x0, 1e4, 1, [], [], cov0, 1);

cov0 = [];

nreps = 4;
niter = [1, 1, 1, 5]*1e4;
for ii = 1:nreps
    [xsto, outsto] = MCMC_adaptive2(obj, x0, niter(ii), 1, cov0, 1);
    inds = find(outsto==max(outsto));
    x0 = xsto(inds(1),:);
    cov0 = cov(xsto);
end




[xsto, outsto] = MCMC_adaptive(obj, x0, niter(ii), 1, [], [], cov0, 1);
save calibration_res;

return;


x2 = xsto2(end,:);
cov0 = cov(xsto2);
[xsto2, outsto2] = MCMC_adaptive(obj, x2, 5e4, 1, [], [], cov0, 1);
fprintf('\n');


nx  = 200;
ix0 = round(size(xsto,1)/2);
dx  = round(size(xsto,1)/2/nx);
xs  = xsto(ix0:dx:end,:);

mk = round(size(xs,1)/24);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ', ii/mk); end 
    [out, aux] = obj(xs(ii,:));
    sim(ii,:) = [aux.incd, aux.mort, aux.p_migrTB, aux.p_migrpopn, aux.p_LTBI];
end
fprintf('\n');

save calibration_res;