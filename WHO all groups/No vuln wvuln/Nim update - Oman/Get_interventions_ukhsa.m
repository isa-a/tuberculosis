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
    nx  = 20;
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

    TPTcov = -log(1-0.25);
    ACFcov = -log(1-0.50);

    % social protection
    % ra1 = r0; pa1 = p0;
    % ra1.TPT = TPTcov * [0 0 0 1];
    % pa1.TPTeff = 0.6;
    % Ma1 = make_model(pa1, ra1, i, s, gps, prm.contmat);    

    % TPT in recent migrants
    ra = r0; pa = p0;
    ra.TPT = TPTcov * [0 1 0 0];
    ra.TPTeff = 0.6;
    Ma = make_model(pa, ra, i, s, gps, prm.contmat);
    
    % TPT at point of entry
    rb = ra; pb = pa;
    pb.migrTPT = 0.25;
    Mb = make_model(pb, rb, i, s, gps, prm.contmat);

    % TPT in longer-term migrants
    % rb1 = rb; pb1 = pb;
    % pb1.TPT = TPTcov * [0 1 1 1];
    % Mb1 = make_model(pb1, rb1, i, s, gps, prm.contmat);

    % TPT in contacts
    % rb2 = rb1; pb2 = pb1;
    % rb2.progression  = rb2.progression * 0.9;
    % rb2.reactivation = rb2.reactivation * 0.9;
    % Mb2 = make_model(pb2, rb2, i, s, gps, prm.contmat);

    % TPT in vulnerables
    rc = rb; pc = pb;
    rc.TPT = TPTcov * [0 1 0 1];
    Mc = make_model(pc, rc, i, s, gps, prm.contmat);

    % ACF in migrants and vulnerables
    % rd = rc; pd = pc;
    % rd.ACF  = ACFcov * [0 1 0 1];
    % rd.ACF2 = rd.ACF;
    % Md = make_model(pd, rd, i, s, gps, prm.contmat);

    % Reducing diagnostic delay
    re = rc; pe = pc;
    re.ACF = ACFcov * [1 0 0 0];
    re.TPT = TPTcov * [0 1 0 1];
    Me = make_model(pe, re, i, s, gps, prm.contmat);

    % New action plan model starting from 2027 with full intervention
    rf = re; pf = pe;
    rf.ACF = -log(1-0.99) * [1 1 1 1];
    rf.TPT = -log(1-0.99) * [1 1 1 1];
    pf.migrTPT = 0.99;
    Mf = make_model(pf, rf, i, s, gps, prm.contmat);

    % Store models: baseline, initial intervention, and action plan
    models = {M0, Me};    
    
    for mi = 1:length(models)
        % Define the ODE system for M0 and Me
        geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, pe, [2024 2027], agg, sel, r0);
        [t, soln] = ode15s(geq, [2022:2041], init, opts);
        
        % Store incidence rates
        sdiff = diff(soln, [], 1);
        incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5;
    end

    % Extract the final state of Me at 2027 as the initial condition for Mf
    idx_2027 = find(t == 2027, 1);
    init_mf = soln(idx_2027, :);  % Final state at 2027 for Mf

    % Run Mf from 2027 to 2041 with scale-up period [2027 2030]
    geq_mf = @(t,in) goveqs_scaleup(t, in, i, s, M0, Mf, p0, pe, [2027 2030], agg, sel, r0);
    [t_mf, soln_mf] = ode15s(geq_mf, [2027:2041], init_mf, opts);

    % Concatenate Me (up to 2026) with Mf (2027 onward) for final incsto
    sdiff_me = sdiff(1:idx_2027-1, i.aux.inc(1));  % Results from Me up to 2026
    sdiff_mf = diff(soln_mf(:, i.aux.inc(1)));     % Results from Mf from 2027 onward
    incsto(:, ii, 3) = [sdiff_me; sdiff_mf] * 1e5;  % Concatenate and scale to 100,000
end
fprintf('\n');

% Define years vector based on the incsto structure
years = 2022:(2022 + size(incsto, 1) - 1);

% Calculate central estimates and percentiles
central_estimate = mean(incsto, 2);              % Mean across samples
lower_bound = prctile(incsto, 2.5, 2);           % 2.5th percentile
upper_bound = prctile(incsto, 97.5, 2);          % 97.5th percentile

% Plotting
ff = figure('Position', [577, 190, 1029, 732]); 
hold on;

% Define colors and legend entries
colors = lines(3);
legendEntries = {'Baseline', 'UKHSA Action Plan', 'New Action Plan'};

% Plot each model with shaded areas
for mi = 1:3
    % Extract the time series for the current model
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    lower_bound_model = squeeze(lower_bound(:, 1, mi));
    upper_bound_model = squeeze(upper_bound(:, 1, mi));

    % Plot shaded area for 2.5th and 97.5th percentiles
    fill([years fliplr(years)], [lower_bound_model' fliplr(upper_bound_model')], ...
        colors(mi, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
    
    % Plot central estimate curve
    plot(years, central_estimate_model, 'LineWidth', 2, 'Color', colors(mi, :)); 
end

% Set axis labels and legend
xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
xlim([years(1) years(end)]);
ylim([0, 8]);
legend(legendEntries, 'FontWeight', 'bold', 'FontSize', 12);

hold off;

