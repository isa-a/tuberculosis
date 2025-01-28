clear all; load Model_setup;

obj = @(x) get_objective(x, prm, ref0, sel0, agg0, gps0, lhd);

load ../../AFR/v0/xx;

obj(xx)

return;

% --- Generate initial guesses

nsam = 1e2;
nx   = 12; %9;

lo = repmat(prm.bounds(1,1:nx),nsam,1);
hi = repmat(prm.bounds(2,1:nx),nsam,1);
sams = lo + (hi-lo).*lhsdesign(nsam,nx);

out = zeros(1,nsam);
inc = zeros(1,nsam);
mk  = round(nsam/25);
for ii = 1:nsam
    if mod(ii,mk)==0; fprintf('%0.5g ',ii/mk); end
    [out(ii), aux, msg(ii)] = obj(sams(ii,:));
    if ~isempty(aux)
        inc(ii) = aux.inc_all;
    end
end
fprintf('\n');

% Order them to get the best-fitting ones
tmp1 = [out;1:nsam]';
tmp2 = tmp1(~isnan(out),:);
mat  = sortrows(tmp2,-1);

ord1  = mat(:,2);
ordx1 = sams(ord1,:);
% ordi = inc(ord1);
% ordo = out(ord1);

% --- Optimise best ones --------------------------------------------------
x0 = ordx1(1:5,:); % ordx1(1:10,:);
% x0 = ordx1(1:2,:);
for ii = 1:size(x0,1)
    [x1(ii,:), fval(ii)] = fminsearch(@(x)-obj(x),x0(ii,:),optimset('PlotFcns',@optimplotfval));    
end

% Order them
tmp1  = [fval; 1:size(x0,1)]';
mat   = sortrows(tmp1,1);
ord2  = mat(:,2);
ordx2 = x1(ord2,:);
% ordo2 = fval(ord2);

load cov0;

% --- Do the sampling -----------------------------------------------------
% xsto = [];
% niter = 1;
% for ii = 1:niter
%     %[xsto(:,:,ii), outsto(ii,:)] = MCMC_adaptive(obj, ordx2(ii,:), 5e4, 1, [], [], cov0, true);
%     [xsto(:,:,ii), outsto(ii,:)] = MCMC_adaptive2(obj, ordx2(ii,:), 5e4, 1, cov0, true);    
% end

x0 = ordx2(ii,:);

niter = 2;
for ii = 1:niter
    %[xsto(:,:,ii), outsto(ii,:)] = MCMC_adaptive(obj, ordx2(ii,:), 5e4, 1, [], [], cov0, true);
    [xsto, outsto] = MCMC_adaptive2(obj, x0, 5e4, 1, cov0, true);    

    xstoall(:,:,ii) = xsto;
    outstoall(:,:,ii) = outsto;

    ind = find(outsto==max(outsto));
    x0   = xsto(ind(1),:);
    cov0 = cov(xsto);
end

save model_fits2;


% ind = 1; mat = xsto(:,:,ind); vec = outsto(ind,:); 
% inds = find(vec==max(vec)); xs = mat(inds(1),:);
% xs = ordx(2,:);
% [out,aux] = obj(xs)
