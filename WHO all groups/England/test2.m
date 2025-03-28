% Initialize storage arrays for outputs
num_models = 9;  % Total number of models
num_time_points = length(2022:2041);  % Fixed years for all models
num_realizations = size(xs, 1);

endsolsto = zeros(num_models, size(init, 2));  % Final solutions for all models
incsto = zeros(num_time_points - 1, num_realizations, num_models);  % Incidence storage
mrtsto = zeros(num_time_points - 1, num_realizations, num_models);  % Mortality storage

% Define chunk size for progress display
mk = round(size(xs, 1) / 25);

for ii = 1:size(xs, 1)
    if mod(ii, mk) == 0
        fprintf('%0.5g ', ii / mk); 
    end

    xx = xs(ii, :);
    [out, aux] = obj(xx);
    init = aux.soln(end, :);

    % Define parameters and initial model M0
    [p0, r0] = allocate_parameters(xx, p, r, xi);
    r0.gamma = r0.gamma_2020;
    M0 = make_model(p0, r0, i, s, gps, prm.contmat);

    % Define intervention parameters for each model
    TPTcov = -log(1 - 0.25);
    ACFcov = -log(1 - 0.50);

    % Define all models
    ra = r0; pa = p0;
    ra.TPT = TPTcov * [0 1 0 0];
    ra.TPTeff = 0.81;
    Ma = make_model(pa, ra, i, s, gps, prm.contmat);

%     rb = ra; pb = pa;
%     pb.migrTPT = 0.25;
%     Mb = make_model(pb, rb, i, s, gps, prm.contmat);

    rc = ra; pc = pa;
    rc.TPT = TPTcov * [0 1 0 1];
    Mc = make_model(pc, rc, i, s, gps, prm.contmat);

    re = rc; pe = pc;
    re.ACF = ACFcov * [1 0 0 0];
    Me = make_model(pe, re, i, s, gps, prm.contmat);

    rf1 = re; pf1 = pe;
    rf1.ACF = -log(1 - 0.99) * [1 1 1 1];
    Mf1 = make_model(pf1, rf1, i, s, gps, prm.contmat);

    rf2 = rf1; pf2 = pf1;
    rf2.TPT = -log(1 - 0.5) * [1 1 1 1];
    pf2.migrTPT = 0.8;
    Mf2 = make_model(pf2, rf2, i, s, gps, prm.contmat);

    rf3 = rf2; pf3 = pf2;
    rf3.progression(:, 4) = rf3.progression(:, 4) / 2;
    rf3.reactivation(:, 4) = rf3.reactivation(:, 4) / 2;
    Mf3 = make_model(pf3, rf3, i, s, gps, prm.contmat);

    rg1 = rf3; pg1 = pf3;
    rg1.TPTeff = 0.9;
    Mg1 = make_model(pg1, rg1, i, s, gps, prm.contmat);

    rg2 = rg1; pg2 = pg1;
    rg2.relapse = rg2.relapse / 2;
    Mg2 = make_model(pg2, rg2, i, s, gps, prm.contmat);

    rg3 = rg2; pg3 = pg2;
    rg3.progression(2, :) = rg3.progression(2, :) * 0.6;
    rg3.reactivation(2, :) = rg3.reactivation(2, :) * 0.6;
    Mg3 = make_model(pg3, rg3, i, s, gps, prm.contmat);


    models = {M0, Ma, Me, Mf1, Mf2, Mf3, Mg1, Mg2, Mg3};

    % Run simulations for all models
    for mi = 1:length(models)
        try
            % Run each model as originally defined
            if mi == 4  % Me: 2024 to 2027
                geq_me = @(t, in) goveqs_scaleup(t, in, i, s, M0, Me, p0, pe, [2024 2027], agg, sel, r0);
                [~, soln] = ode15s(geq_me, 2022:2041, init, opts);

            elseif mi >= 5 && mi <= 7  % Mf1, Mf2, Mf3: 2027 to 2030
                geq_me = @(t, in) goveqs_scaleup(t, in, i, s, M0, Me, p0, pe, [2024 2027], agg, sel, r0);
                [~, soln1] = ode15s(geq_me, 2022:2027, init, opts);

                geq_mf = @(t, in) goveqs_scaleup(t, in, i, s, Me, models{mi}, pe, pf1, [2027 2030], agg, sel, r0);
                [~, soln2] = ode15s(geq_mf, 2027:2041, soln1(end, :), opts);

                soln = [soln1; soln2(2:end, :)];

            elseif mi >= 8 && mi <= 9  % Mg2, Mg3: 2030 to 2033
                geq_me = @(t, in) goveqs_scaleup(t, in, i, s, M0, Me, p0, pe, [2024 2027], agg, sel, r0);
                [~, soln1] = ode15s(geq_me, 2022:2027, init, opts);

                geq_mf = @(t, in) goveqs_scaleup(t, in, i, s, Me, Mf3, pe, pf3, [2027 2030], agg, sel, r0);
                [~, soln2] = ode15s(geq_mf, 2027:2030, soln1(end, :), opts);

                geq_mg = @(t, in) goveqs_scaleup(t, in, i, s, Mf3, models{mi}, pf3, pg2, [2030 2033], agg, sel, r0);
                [~, soln3] = ode15s(geq_mg, 2030:2041, soln2(end, :), opts);

                soln = [soln1; soln2(2:end, :); soln3(2:end, :)];

            else
                % Standard simulation
                geq = @(t, in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, pe, [2024 2027], agg, sel, r0);
                [~, soln] = ode15s(geq, 2022:2041, init, opts);
            end

            % Store results
            endsolsto(mi, :) = soln(end, :);
            sdiff = diff(soln, [], 1);
            incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5;
            mrtsto(:, ii, mi) = sdiff(:, i.aux.mort) * 1e5;

        catch ME
            fprintf('Error in model %d: %s\n', mi, ME.message);
        end
    end
end
fprintf('\n');


% Plotting logic remains the same as the original code


% Plotting code remains unchanged except `Mg1` is excluded from the legend setup

% Define the years and calculate central estimates and percentiles
years = 2022:(2022 + size(incsto, 1) - 1);
central_estimate = mean(incsto, 2);              
lower_bound = prctile(incsto, 2.5, 2);          
upper_bound = prctile(incsto, 97.5, 2);          

% Create figure
ff = figure('Position', [577, 190, 1029, 732]); 
hold on;

% Set colors for model groups
baseline_color = [0 0 1];                        % Blue for Baseline (M0)
ma_me_color = [1 0, 0];                           % Red for Ma and Me
mf1_mf2_mf3_color = [1 0.5, 0];                   % Orange for Mf1, Mf2, Mf3
mg1_mg2_mg3_color = [0.4940, 0.1840, 0.5560];     % MATLAB purple for Mg1, Mg2, Mg3

% Define the legend entries
legendEntries = {'Baseline', 'Ma', 'Me', 'Mf1', 'Mf2', 'Mf3', 'Mg1', 'Mg2', 'Mg3'};

% Initialize plot handle array
plotHandles = gobjects(length(models), 1);

% Iterate over models
for mi = 1:length(models)
    % Assign color based on model group
    if mi == 1  % Baseline (M0)
        plot_color = baseline_color;
    elseif mi == 2 || mi == 3  % Ma or Me
        plot_color = ma_me_color;
    elseif mi >= 4 && mi <= 6  % Mf1, Mf2, Mf3
        plot_color = mf1_mf2_mf3_color;
    elseif mi >= 7 && mi <= 9  % Mg1, Mg2, Mg3
        plot_color = mg1_mg2_mg3_color;
    end

    % Assign line style
    if mi == 1 || mi == 3 || mi == 6 || mi == 9  % Baseline, Me, Mf3, Mg3
        line_style = '-';  % Solid
    else
        line_style = '--';  % Dashed
    end

    % Extract data for the current model
    central_estimate_model = central_estimate(:, mi);
    lower_bound_model = lower_bound(:, mi);
    upper_bound_model = upper_bound(:, mi);

    % Add shading for Baseline, Me, Mf3, Mg3
    if mi == 1 || mi == 3 || mi == 6 || mi == 9  % Baseline, Me, Mf3, Mg3
        fill([years fliplr(years)], ...
             [lower_bound_model' fliplr(upper_bound_model')], ...
             plot_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    % Plot the central estimate
    plotHandles(mi) = plot(years, central_estimate_model, 'LineWidth', 2, ...
                           'Color', plot_color, 'LineStyle', line_style);
end

% Add legend
legend(plotHandles, legendEntries, 'FontWeight', 'bold', 'FontSize', 12, 'Location', 'best');

% Add labels and axis limits
xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
xlim([years(1) years(end)]);
ylim([0, 10]);

% Add horizontal reference lines
yline(1, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
yline(0.1, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');

hold off;

