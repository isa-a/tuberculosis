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

    % TPT in recent migrants
    ra = r0; pa = p0;
    ra.TPT = TPTcov * [0 1 0 0];
    ra.TPTeff = 0.6;
    Ma = make_model(pa, ra, i, s, gps, prm.contmat);
    
    % TPT at point of entry
    rb = ra; pb = pa;
    pb.migrTPT = 0.25;
    Mb = make_model(pb, rb, i, s, gps, prm.contmat);

    % TPT in vulnerables
    rc = rb; pc = pb;
    rc.TPT = TPTcov * [0 1 0 1];
    Mc = make_model(pc, rc, i, s, gps, prm.contmat);

    % Reducing diagnostic delay
    re = rc; pe = pc;
    re.ACF = ACFcov * [1 0 0 0];
    re.TPT = TPTcov * [0 1 0 1];
    Me = make_model(pe, re, i, s, gps, prm.contmat);

    % Full intervention plan (`Mf`)
    rf = re; pf = pe;
    rf.ACF = -log(1-0.99) * [1 1 1 1];
    rf.TPT = -log(1-0.99) * [1 1 1 1];
    pf.migrTPT = 0.99;
    Mf = make_model(pf, rf, i, s, gps, prm.contmat);

    % New action plan (`Mg`) starting from 2030
    rg = re; pg = pe;
    rg.ACF = -log(1-0.99) * [1 1 1 1];
    rg.TPT = -log(1-0.99) * [1 1 1 1];
    pg.migrTPT = 0.99;
    rg.TPTeff = 0.8;
    rg.relapse = rg.relapse/2;
    rg.progression = rg.progression/4;
    rg.reactivation = rg.reactivation/4; 
    Mg = make_model(pg, rg, i, s, gps, prm.contmat);

    % Store models: baseline, initial intervention, action plan, and new action plan
    models = {M0, Me, Mf, Mg};    
    
    for mi = 1:length(models)
        if mi == length(models)
            % Run Mf up to 2030, then continue with Mg from 2030 onward
            geq_mf = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi-1}, p0, pf, [2024 2030], agg, sel, r0);
            [t1, soln1] = ode15s(geq_mf, 2022:2030, init, opts);

            % Final solution for Mf at 2030
            init_final = soln1(end,:);

            % Mg from 2030 onwards
            geq_mg = @(t,in) goveqs_scaleup(t, in, i, s, models{mi-1}, models{mi}, p0, pg, [2030 2033], agg, sel, r0);
            [t2, soln2] = ode15s(geq_mg, 2030:2041, init_final, opts);

            % Combine the solutions
            t = [t1; t2(2:end)];
            soln = [soln1; soln2(2:end,:)];
        elseif mi == length(models)-1
            % For Mf, run Me up to 2027, then continue with Mf from 2027 onward
            geq_me = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi-1}, p0, pe, [2024 2027], agg, sel, r0);
            [t1, soln1] = ode15s(geq_me, 2022:2027, init, opts);

            % Final solution for Me at 2027
            init_final = soln1(end,:);

            % Mf from 2027 onwards
            geq_mf = @(t,in) goveqs_scaleup(t, in, i, s, models{mi-1}, models{mi}, p0, pf, [2027 2030], agg, sel, r0);
            [t2, soln2] = ode15s(geq_mf, 2027:2041, init_final, opts);

            % Combine the solutions
            t = [t1; t2(2:end,:)];
            soln = [soln1; soln2(2:end,:)];
        else
            % Every other scenario remains the same
            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, pa, [2024 2029], agg, sel, r0);
            [t, soln] = ode15s(geq, 2022:2041, init, opts);
        end

        endsolsto(mi,:) = soln(end,:);

        sdiff = diff(soln, [], 1);
        incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5;
        mrtsto(:, ii, mi) = sdiff(:, i.aux.mort) * 1e5;
    end
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
colors = lines(length(models));
legendEntries = {'Baseline', 'UKHSA Action Plan', 'Full Intervention Plan', 'New Action Plan'};

% Initialize an array to store plot handles for the legend
plotHandles = gobjects(length(models),1);

% Plot each model with shaded areas
for mi = 1:length(models)
    % Extract the time series for the current model
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    lower_bound_model = squeeze(lower_bound(:, 1, mi));
    upper_bound_model = squeeze(upper_bound(:, 1, mi));
    
    if mi == length(models)
        % For Mg, use Mf's estimates up to 2030
        idx_2030 = find(years == 2030);
        
        % Use Mf's estimates up to 2030
        central_estimate_model(1:idx_2030) = central_estimate(1:idx_2030, 1, mi-1);
        lower_bound_model(1:idx_2030) = lower_bound(1:idx_2030, 1, mi-1);
        upper_bound_model(1:idx_2030) = upper_bound(1:idx_2030, 1, mi-1);
    end

    % Plot shaded area for 2.5th and 97.5th percentiles
    fillHandle = fill([years fliplr(years)], [lower_bound_model' fliplr(upper_bound_model')], ...
        colors(mi, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
    
    % Plot central estimate curve
    plotHandle = plot(years, central_estimate_model, 'LineWidth', 2, 'Color', colors(mi, :)); 
    
    % Store the plot handle for the legend
    plotHandles(mi) = plotHandle;
    
    % Update the legend to include only the plotted models so far
    legend(plotHandles(1:mi), legendEntries(1:mi), 'FontWeight', 'bold', 'FontSize', 12);
    
    % Set axis labels and limits (only need to set once, but harmless to set multiple times)
    xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
    xlim([years(1) years(end)]);
    ylim([0, 9]);
    
    % Pause execution and wait for the user to press a key
    disp('Press the spacebar to display the next curve...');
    pause;  % Waits for the user to press a key
end

hold off;