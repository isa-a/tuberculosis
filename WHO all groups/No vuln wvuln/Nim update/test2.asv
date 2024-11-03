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

    % Reducing diagnostic delay (Me)
    re = rc; pe = pc;  % Start from baseline parameters
    re.ACF = ACFcov * [1 0 0 0];
    re.TPT = TPTcov * [0 1 0 1];
    pe.migrTPT = 0.25;  % Include TPT at point of entry for migrants
    Me = make_model(pe, re, i, s, gps, prm.contmat);

    % 2027
    rf = re; pf = pe;
    rf.ACF = -log(1-0.99) * [1 1 1 1];
    rf.TPT = -log(1-0.99) * [1 1 1 1];
    pf.migrTPT = 0.8;
    rf.progression(:,4) = rf.progression(:,4)/2;
    rf.reactivation(:,4) = rf.reactivation(:,4)/2;
    Mf = make_model(pf, rf, i, s, gps, prm.contmat);

    % 2030
    rg = rf; pg = pf;
    rg.ACF = -log(1-0.99) * [1 1 1 1];
    rg.TPT = -log(1-0.99) * [1 1 1 1];
    pg.migrTPT = 0.8;
    rg.TPTeff = 0.8;
    rg.relapse = rg.relapse/2;
    rg.progression(2,:) = rg.progression(2,:)*0.6;
    rg.reactivation(2,:) = rg.reactivation(2,:)*0.6; 
    Mg = make_model(pg, rg, i, s, gps, prm.contmat);

    % Include Ma, Mb, Mc in models array
    models = {M0, Ma, Mb, Mc, Me, Mf, Mg};    
    
    for mi = 1:length(models)
        if mi == length(models)
            % mi == 7, Mg
            % Run Mf to 2030, then Mg
            % First, run Me up to 2027
            geq_me = @(t,in) goveqs_scaleup(t, in, i, s, M0, Me, p0, pe, [2024 2027], agg, sel, r0);
            [t1, soln1] = ode15s(geq_me, 2022:2027, init, opts);

            % Then, run Mf from 2027 to 2030
            geq_mf = @(t,in) goveqs_scaleup(t, in, i, s, Me, Mf, pe, pf, [2027 2030], agg, sel, r0);
            [t2, soln2] = ode15s(geq_mf, 2027:2030, soln1(end,:), opts);

            % Finally, run Mg from 2030 onwards
            geq_mg = @(t,in) goveqs_scaleup(t, in, i, s, Mf, Mg, pf, pg, [2030 2033], agg, sel, r0);
            [t3, soln3] = ode15s(geq_mg, 2030:2041, soln2(end,:), opts);

            % Combine
            t = [t1; t2(2:end); t3(2:end)];
            soln = [soln1; soln2(2:end,:); soln3(2:end,:)];
        elseif mi == length(models)-1
            % mi == 6, Mf
            % For Mf, run Me up to 2027, then Mf 
            geq_me = @(t,in) goveqs_scaleup(t, in, i, s, M0, Me, p0, pe, [2024 2027], agg, sel, r0);
            [t1, soln1] = ode15s(geq_me, 2022:2027, init, opts);

            % Mf from 2027
            geq_mf = @(t,in) goveqs_scaleup(t, in, i, s, Me, Mf, pe, pf, [2027 2030], agg, sel, r0);
            [t2, soln2] = ode15s(geq_mf, 2027:2041, soln1(end,:), opts);

            % Combine
            t = [t1; t2(2:end)];
            soln = [soln1; soln2(2:end,:)];
        elseif mi == 5
            % mi == 5, Me
            % Simulate Me starting from M0
            geq_me = @(t,in) goveqs_scaleup(t, in, i, s, M0, Me, p0, pe, [2024 2027], agg, sel, r0);
            [t, soln] = ode15s(geq_me, 2022:2041, init, opts);
        else
            % All others (including Ma, Mb, Mc)
            % Simulate from M0 to the intervention model
            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, pe, [2024 2027], agg, sel, r0);
            [t, soln] = ode15s(geq, 2022:2041, init, opts);
        end

        % Store results
        endsolsto(mi,:) = soln(end,:);

        sdiff = diff(soln, [], 1);
        incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5;
        mrtsto(:, ii, mi) = sdiff(:, i.aux.mort) * 1e5;
    end
end
fprintf('\n');

% Prepare for plotting
years = 2022:(2022 + size(incsto, 1) - 1);

% Percentiles
central_estimate = mean(incsto, 2);              
lower_bound = prctile(incsto, 2.5, 2);          
upper_bound = prctile(incsto, 97.5, 2);          

% Plot
ff = figure('Position', [577, 190, 1029, 732]); 
hold on;

% Generate colors
colors = lines(length(models));

% Set Me and Ma, Mb, Mc to red
colors(5,:) = [1 0 0];  % Me
colors(2,:) = colors(5,:);  % Ma
colors(3,:) = colors(5,:);  % Mb
colors(4,:) = colors(5,:);  % Mc

% Define legend entries
legendEntries = {'Baseline', 'Recent migrants TPT', 'Migrants TPT at entry', 'Vulnerable TPT', 'UKHSA Action Plan', 'Max scaleup of new tools', 'New TPT and treatments'};

% Initialize plotHandles
plotHandles = gobjects(length(models),1);

for mi = 1:length(models)
    % Extract data for model mi
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    lower_bound_model = squeeze(lower_bound(:, 1, mi));
    upper_bound_model = squeeze(upper_bound(:, 1, mi));

    % Adjust data for Mf to start from 2027
    if mi == length(models)-1  % mi == 6, Mf
        idx_2027 = find(years == 2027);

        central_estimate_model(1:idx_2027-1) = central_estimate(1:idx_2027-1, 1, mi-1);
        lower_bound_model(1:idx_2027-1) = lower_bound(1:idx_2027-1, 1, mi-1);
        upper_bound_model(1:idx_2027-1) = upper_bound(1:idx_2027-1, 1, mi-1);
    end

    % Adjust data for Mg to start from 2030
    if mi == length(models)  % mi == 7, Mg
        idx_2030 = find(years == 2030);

        central_estimate_model(1:idx_2030-1) = central_estimate(1:idx_2030-1, 1, mi-1);
        lower_bound_model(1:idx_2030-1) = lower_bound(1:idx_2030-1, 1, mi-1);
        upper_bound_model(1:idx_2030-1) = upper_bound(1:idx_2030-1, 1, mi-1);
    end

    % Set plot color and line style
    if mi == 5  % Me
        plot_color = colors(mi, :);
        line_style = '-';
        plot_shading = true;
    elseif mi >=2 && mi <=4  % Ma, Mb, Mc
        plot_color = colors(mi, :);  % Red color
        line_style = '--';  % Dashed line
        plot_shading = false;  % No shading
    elseif mi == 1  % Baseline
        plot_color = colors(mi, :);
        line_style = '-';
        plot_shading = true;  % Shading for baseline
    else
        plot_color = colors(mi, :);
        line_style = '-';
        plot_shading = true;  % Shading for Mf and Mg
    end

    % Plot fill area (uncertainty) if applicable
    if plot_shading
        % Plot fill area (uncertainty)
        fillHandle = fill([years fliplr(years)], [lower_bound_model' fliplr(upper_bound_model')], ...
            plot_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    % Plot central estimate
    plotHandle = plot(years, central_estimate_model, 'LineWidth', 2, 'Color', plot_color, 'LineStyle', line_style);

    % Store plot handle for legend
    plotHandles(mi) = plotHandle;

    % Update legend
    legend(plotHandles(1:mi), legendEntries(1:mi), 'FontWeight', 'bold', 'FontSize', 12);

    % Set labels and limits
    xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
    xlim([years(1) years(end)]);
    ylim([0, 9]);
    yline(1,'k--','LineWidth', 2, 'HandleVisibility', 'off');
    yline(0.1,'k--','LineWidth', 2, 'HandleVisibility', 'off');

    pause;  % Waits for the user to press a key
end

hold off;
