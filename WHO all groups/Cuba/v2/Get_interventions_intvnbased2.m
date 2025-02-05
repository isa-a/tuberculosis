% Determine model iteration checkpoints
mk = round(size(xs,1)/25);

% Initialize result storage for each simulation run
n_years = 2041 - 2022 + 1;  % Duration from 2022 to 2041
num_models = 9;  % Total number of models: M0, Ma, Me, Mf1, Mf2, Mf3, Mg1, Mg2, Mg3
incsto = nan(n_years, size(xs, 1), num_models);  % Initialize with NaN
mrtsto = nan(n_years, size(xs, 1), num_models);  % Initialize with NaN
endsolsto = nan(num_models, size(init, 2), size(xs, 1));  % Initialize endsolsto for each model and simulation

% Loop through simulations
for ii = 1:size(xs,1)
    
    if mod(ii, mk) == 0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out, aux] = obj(xx);
    init = aux.soln(end,:);

    % Define parameters and initial model M0
    [p0, r0] = allocate_parameters(xx, p, r, prm, xi);
    r0.gamma = r0.gamma_2020;
    M0 = make_model(p0, r0, i, s, gps, prm.contmat);

    % Define intervention parameters for each model
    TPTcov = -log(1 - 0.25);
    ACFcov = -log(1 - 0.50);

    % Define all models with appropriate parameters
    ra = r0; pa = p0;
    ra.TPT = TPTcov * [0 1 0 0];
    ra.TPTeff = 0.6;
    Ma = make_model(pa, ra, i, s, gps, prm.contmat);

    re = ra; pe = pa; 
    re.ACF = -log(1 - 0.99) * [1 1 1 1];
    Me = make_model(pe, re, i, s, gps, prm.contmat);

    rf1 = re; pf1 = pe;
    rf1.ACF = -log(1 - 0.99) * [1 1 1 1];
    rf1.TPT = -log(1 - 0.5) * [1 1 1 1];
    Mf1 = make_model(pf1, rf1, i, s, gps, prm.contmat);

    rf2 = rf1; pf2 = pf1;
    pf2.migrTPT = 0.8;
    Mf2 = make_model(pf2, rf2, i, s, gps, prm.contmat);

    rf3 = rf2; pf3 = pf2;
    rf3.progression(:, 4) = rf3.progression(:, 4) / 2;
    rf3.reactivation(:, 4) = rf3.reactivation(:, 4) / 2;
    Mf3 = make_model(pf3, rf3, i, s, gps, prm.contmat);

    rg1 = rf3; pg1 = pf3;
    rg1.TPTeff = 0.8;
    Mg1 = make_model(pg1, rg1, i, s, gps, prm.contmat);

    rg2 = rg1; pg2 = pg1;
    rg2.relapse = rg2.relapse / 2;
    Mg2 = make_model(pg2, rg2, i, s, gps, prm.contmat);

    rg3 = rg2; pg3 = pg2;
    rg3.progression(2, :) = rg3.progression(2, :) * 0.6;
    rg3.reactivation(2, :) = rg3.reactivation(2, :) * 0.6;
    Mg3 = make_model(pg3, rg3, i, s, gps, prm.contmat);

    % Collect all models in array
    models = {M0, Ma, Me, Mf1, Mf2, Mf3, Mg1, Mg2, Mg3};    
    
    % Run simulations with updated scale-up periods and bifurcation points
    for mi = 1:length(models)
        if mi == 2  % Ma: 2024 to 2027
            geq_ma = @(t, in) goveqs_scaleup(t, in, i, s, M0, Ma, p0, pa, [2024 2027], agg, sel, r0);
            [t, soln] = ode15s(geq_ma, 2022:2041, init, opts);

        elseif mi == 3  % Me: 2024 to 2027
            geq_me = @(t, in) goveqs_scaleup(t, in, i, s, M0, Me, p0, pe, [2024 2027], agg, sel, r0);
            [t, soln] = ode15s(geq_me, 2022:2041, init, opts);

        elseif mi >= 4 && mi <= 6  % Mf1, Mf2, Mf3: 2027 to 2030, bifurcating from Me
            geq_me = @(t, in) goveqs_scaleup(t, in, i, s, M0, Me, p0, pe, [2024 2027], agg, sel, r0);
            [t1, soln1] = ode15s(geq_me, 2022:2027, init, opts);
            
            % Select appropriate parameter set for each Mf model
            if mi == 4
                geq_mf = @(t, in) goveqs_scaleup(t, in, i, s, Me, Mf1, pe, pf1, [2027 2030], agg, sel, rf1);
            elseif mi == 5
                geq_mf = @(t, in) goveqs_scaleup(t, in, i, s, Me, Mf2, pe, pf2, [2027 2030], agg, sel, rf2);
            elseif mi == 6
                geq_mf = @(t, in) goveqs_scaleup(t, in, i, s, Me, Mf3, pe, pf3, [2027 2030], agg, sel, rf3);
            end
            [t2, soln2] = ode15s(geq_mf, 2027:2041, soln1(end, :), opts);

            t = [t1; t2(2:end)];
            soln = [soln1; soln2(2:end, :)];

        elseif mi >= 7 && mi <= 9  % Mg1, Mg2, Mg3: 2030 to 2033, bifurcating from Mf3
            geq_me = @(t, in) goveqs_scaleup(t, in, i, s, M0, Me, p0, pe, [2024 2027], agg, sel, r0);
            [t1, soln1] = ode15s(geq_me, 2022:2027, init, opts);
            
            geq_mf = @(t, in) goveqs_scaleup(t, in, i, s, Me, Mf3, pe, pf3, [2027 2030], agg, sel, rf3);
            [t2, soln2] = ode15s(geq_mf, 2027:2030, soln1(end, :), opts);

            % Select appropriate parameter set for each Mg model
            if mi == 7
                geq_mg = @(t, in) goveqs_scaleup(t, in, i, s, Mf3, Mg1, pf3, pg1, [2030 2033], agg, sel, rg1);
            elseif mi == 8
                geq_mg = @(t, in) goveqs_scaleup(t, in, i, s, Mf3, Mg2, pf3, pg2, [2030 2033], agg, sel, rg2);
            elseif mi == 9
                geq_mg = @(t, in) goveqs_scaleup(t, in, i, s, Mf3, Mg3, pf3, pg3, [2030 2033], agg, sel, rg3);
            end
            [t3, soln3] = ode15s(geq_mg, 2030:2041, soln2(end, :), opts);

            t = [t1; t2(2:end); t3(2:end)];
            soln = [soln1; soln2(2:end, :); soln3(2:end, :)];
        else
            % Baseline model
            geq = @(t, in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, pe, [2024 2027], agg, sel, r0);
            [t, soln] = ode15s(geq, 2022:2041, init, opts);
        end

        endsolsto(mi, :, ii) = soln(end, :);  % Store solutions with () indexing
        sdiff = diff(soln, [], 1);
        incsto(1:size(sdiff, 1), ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5;
        mrtsto(1:size(sdiff, 1), ii, mi) = sdiff(:, i.aux.mort) * 1e5;
    end
end
fprintf('\n');

% Plotting
% Define the years and calculate central estimates and percentiles
% Plotting
% Define the years and calculate central estimates and percentiles
years = 2022:(2022 + size(incsto, 1) - 1);
central_estimate = mean(incsto, 2, 'omitnan');              
lower_bound = prctile(incsto, 2.5, 2);          
upper_bound = prctile(incsto, 97.5, 2);          

figure('Position', [577, 190, 1029, 732]); 
hold on;

% Plot styles
baseline_color = [0 0 1];  
ma_me_color = [1 0 0];  
mf1_mf2_mf3_color = [1 0.5 0];  
mg1_mg2_mg3_color = [0.5 0 0.5];
solid_models = [1, 3, 6, 9];
legendEntries = {'Baseline', 'maT', 'me', 'Mf1', 'Mf2', 'Mf3', 'Mg1', 'Mg2', 'Mg3'};
plotHandles = gobjects(length(models),1);

% Plot each model
for mi = 1:length(models)
    central_estimate_model = squeeze(central_estimate(:,1,mi));
    lower_bound_model = squeeze(lower_bound(:,1,mi));
    upper_bound_model = squeeze(upper_bound(:,1,mi));
    
    % Determine plot color based on model
    if mi == 1
        plot_color = baseline_color; 
    elseif ismember(mi, [2, 3])
        plot_color = ma_me_color;
    elseif ismember(mi, [4, 5, 6])
        plot_color = mf1_mf2_mf3_color;
    elseif ismember(mi, [7, 8, 9])
        plot_color = mg1_mg2_mg3_color;
    end
    
    % If the model is one of the specified models, add shading
    if ismember(mi, [1, 3, 6, 9])
        % Create coordinates for the shaded area
        x_shade = [years, fliplr(years)];
        y_shade = [upper_bound_model', fliplr(lower_bound_model')];
        
        % Plot the shaded area
        fill_handle = fill(x_shade, y_shade, plot_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        uistack(fill_handle, 'bottom');  % Ensure shading is behind the lines
    end
    
    % Determine line style
    line_style = '-'; 
    if ~ismember(mi, solid_models)
        line_style = '--'; 
    end
    
    % Plot the central estimate line
    plotHandle = plot(years, central_estimate_model, 'LineWidth', 2, 'Color', plot_color, 'LineStyle', line_style);
    plotHandles(mi) = plotHandle;
end

% Finalize the plot
legend(plotHandles, legendEntries, 'FontWeight', 'bold', 'FontSize', 12);
xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
xlim([years(1) 2040]);
ylim([0, 9]);
yline(1, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
yline(0.1, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');

hold off;

