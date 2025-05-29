clear all; load optim_res1000.mat;

% rin_vec(end-2:end) = rin_vec(end-3);
obj = @(x) get_objective3(x, ref, prm, gps, prm.contmat, rin_vec, lhd);

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

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
    incd(ii) = aux.incd(end);

    [p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);
    p0.prev_in_migr = 0;
    r0.gamma        = r0.gamma_2015;
    M0              = make_model(p0,r0,i,s,gps,prm0.contmat);

    % >2015: scaleup of TPT
    p1 = p0; r1 = r0; prm1 = prm0;
    r1.TPT   = [0 r0.TPT2020rec 0];
    r1.gamma = r0.gamma_2015;
    M1       = make_model(p1,r1,i,s,gps,prm1.contmat);

    % >2010: increase in case-finding
    p2 = p0; r2 = r0; prm2 = prm0;
    p2.prev_in_migr = 0;
    r2.gamma        = r0.gamma_2020;
    M2              = make_model(p2,r2,i,s,gps,prm2.contmat);

    % UKHSA action plan (from now)
    TPTcov = -log(1-0.25);
    ACFcov = -log(1-0.50); 
    ra = r0; pa = p0; prma = prm0;
    ra.TPT = TPTcov * [0 1 0 0];
    ra.ACF = ACFcov * [1 0 0 0];
    Ma = make_model(pa, ra, i, s, gps, prma.contmat);

    % Expanded deployment of current tools (2027 – 2030)
    rb = ra; pb = pa; prmb = prma;
    rb.ACF = -log(1-0.99) * [1 1 1 1];
    rb.TPT = -log(1-0.5)  * [1 1 1 1];
    pb.migrTPT = 0.8;
    Mb         = make_model(pb, rb, i, s, gps, prmb.contmat);

    % same as in obj
%     geq = @(t,in) goveqs_basis3(t, in, i, s, M0, rin_vec, agg, sel, r0, p0, false);
%     [t, soln] = ode15s(geq, 2022:2041, init, opts);    
    
    geq1 = @(t,in) goveqs_scaleup2D(t, in, M0, M1, M2, rin_vec, [2015 2020; 2010 2011], i, s, p2, r2, prm0, sel, agg, false);
    [t1, soln1] = ode15s(geq1, 2022:2041, init, opts);
    sdiff = diff(soln1, [], 1);
    pops  = sum(soln1(:,1:i.nstates),2);
    incsto(:, ii) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
end
fprintf('\n');

mat = prctile(incsto,[2.5,50,97.5],2);
prctile(incd,[2.5,50,97.5]);

figure; plot(mat);

return;

% years vector based on incsto
years = 2022:(2022 + size(incsto, 1) - 1);

central_estimate = mean(incsto, 2);             
lowerbound = prctile(incsto, 2.5, 2);          
upperbound = prctile(incsto, 97.5, 2);        

ff = figure('Position', [577, 190, 1029, 732]); 
hold on;

colors = lines(3);
legendEntries = {'Baseline', 'UKHSA Action Plan', 'New Action Plan'};

for mi = 1:3
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    lower = squeeze(lowerbound(:, 1, mi));
    upper = squeeze(upperbound(:, 1, mi));
    fill([years fliplr(years)], [lower' fliplr(upper')], colors(mi, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(years, central_estimate_model, 'LineWidth', 2, 'Color', colors(mi, :));
end

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
xlim([years(1) years(end)]);
ylim([0, 11]);
legend(legendEntries, 'FontWeight', 'bold', 'FontSize', 12);

hold off;
