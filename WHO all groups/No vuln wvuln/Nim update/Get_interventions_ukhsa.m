clear all; load optim_res_MAIN.mat; load Model_setup;

obj = @(x) get_objective2(x, ref, prm, gps, prm.contmat, lhd);

set1 = {'ds'};
set2 = {'dom','migr','vuln'};
set3 = {'L','P','R','T'};
[inci, incs, incd, lim] = get_addresses({set3, set2, set1}, [], [], [], 0);
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

midpt = false; 
if midpt
    xs = x0sto(2,:);
else
    ix0 = size(xsto,1)/2;
    nx  = 200;
    dx  = round(ix0/nx);
    xs  = xsto(ix0:dx:end,:);
end

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);
    
    init = aux.soln(end,:);

    [p0,r0] = allocate_parameters(xx,p,r,xi);
    r0.gamma = r0.gamma_2020;
    M0 = make_model(p0,r0,i,s,gps,prm.contmat);

    TPTcov = -log(1-0.99);
    ACFcov = -log(1-0.99);

    ra1 = r0; pa1 = p0;
    ra1.TPT = TPTcov*[0 0 0 1];
    pa1.TPTeff = 0.8;
    Ma1 = make_model(pa1,ra1,i,s,gps,prm.contmat);    

    ra = ra1; pa = pa1;
    ra.TPT = TPTcov*[0 1 0 0];
    Ma = make_model(pa,ra,i,s,gps,prm.contmat);
    
    rb = ra; pb = pa;
    pb.migrTPT = 1;
    Mb = make_model(pb,rb,i,s,gps,prm.contmat);

    rb2 = rb; pb2 = pb;
    rb2.progression  = rb2.progression*0.9;
    rb2.reactivation = rb2.reactivation*0.9;
    Mb2 = make_model(pb2,rb2,i,s,gps,prm.contmat);

    rc = rb2; pc = pb2;
    rc.TPT = TPTcov*[0 1 1 1];
    Mc = make_model(pc,rc,i,s,gps,prm.contmat);

    rd = rc; pd = pc;
    rd.ACF  = ACFcov*[0 1 0 1];
    rd.ACF2 = rd.ACF;
    Md = make_model(pd,rd,i,s,gps,prm.contmat);

    re = rd; pe = pd;
    re.ACF  = ACFcov*[1 1 1 1];
    re.ACF2 = rd.ACF;
    re.TPT  = TPTcov*[1 1 1 1];
    Me = make_model(pe,re,i,s,gps,prm.contmat);

    models = {M0, Ma1, Ma, Mb, Mb2, Mc, Md, Me};    
    for mi = 1:length(models)
        
        geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, pa, [2024 2029], agg, sel, r0);
        [t,soln] = ode15s(geq, [2022:2041], init, opts);
        
        sdiff = diff(soln,[],1);
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1))*1e5;
    end
end
fprintf('\n');

years = 2022:(2022 + size(incsto, 1) - 1);

% Compute central estimate and percentiles across samples (2nd dimension)
central_estimate = mean(incsto, 2);               % Mean across samples
lower_bound = prctile(incsto, 2.5, 2);            % 2.5th percentile across samples
upper_bound = prctile(incsto, 97.5, 2);           % 97.5th percentile across samples

% Define a color map for each curve and shading
colors = lines(length(models));  % MATLAB's 'lines' colormap for distinct colors

ff = figure('Position', [577, 190, 1029, 732]); 
lw = 1.5; 
fs = 14;
hold on;

% Plot shaded areas and curves for each scenario
for mi = 1:length(models)
    % Extract the time series for the current model
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    lower_bound_model = squeeze(lower_bound(:, 1, mi));
    upper_bound_model = squeeze(upper_bound(:, 1, mi));

    % Use the color for both the shaded area and the line
    color = colors(mi, :);

    % Plot the shaded area for 2.5th and 97.5th percentiles, excluding from legend
    fill([years fliplr(years)], [lower_bound_model' fliplr(upper_bound_model')], ...
        color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
    
    % Plot the central estimate curve (included in legend)
    plot(years, central_estimate_model, 'LineWidth', 2, 'Color', color); 
end

% Axis labels and settings
xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'FontWeight', 'bold', 'FontSize', 12);
xlim([years(1) years(end)]);
legend(legendEntries, 'FontWeight', 'bold', 'FontSize', 12);

hold off;


