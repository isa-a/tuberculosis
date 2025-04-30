clear all; load mcmc_res.mat;

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
incsto = zeros(20, size(xs,1), 6); % 6 models, 20 years
for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);
      
    init = aux.init;

    [p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);
    p0.prev_in_migr = 0;
    r0.gamma = r0.gamma_2015;
    r0.TPT = [0 r0.TPT2020rec 0];
    M0 = make_model(p0, r0, i, s, gps, prm0.contmat);

    TPTcov = -log(1-0.25);
    ACFcov = -log(1-0.50); 

    % Expanded deployment (2027–2030)
    rb = r0; pb = p0; prmb = prm0;
    rb.ACF = -log(1-0.99) * [1 1 1 1];
    rb.TPT = -log(1-0.5) * [0 1 0 0];
    pb.migrTPT = 0.8;
    Mb = make_model(pb, rb, i, s, gps, prmb.contmat);

    % Mb: ACF domestic only
    rb_ACFdom = r0; pb_ACFdom = p0; prmb_ACFdom = prm0;
    rb_ACFdom.ACF = -log(1-0.99) * [1 0 0 0];
    Mb_ACFdom = make_model(pb_ACFdom, rb_ACFdom, i, s, gps, prmb_ACFdom.contmat);

    % Mb: ACF migrant and domestic
    rb_ACFmig = r0; pb_ACFmig = p0; prmb_ACFmig = prm0;
    rb_ACFmig.ACF = -log(1-0.99) * [1 1 0 0];
    Mb_ACFmig = make_model(pb_ACFmig, rb_ACFmig, i, s, gps, prmb_ACFmig.contmat);

    % Mb: and TPT migrant 
    rb_TPT = r0; pb_TPT = p0; prmb_TPT = prm0;
    rb_TPT.ACF = -log(1-0.99) * [1 1 0 0];
    rb_TPT.TPT = -log(1-0.5) * [0 1 0 0];
    Mb_TPT = make_model(pb_TPT, rb_TPT, i, s, gps, prmb_TPT.contmat);

    % Mb: and migrant entry TPT 
    rb_migTPT = r0; pb_migTPT = p0; prmb_migTPT = prm0;
    rb_migTPT.ACF = -log(1-0.99) * [1 1 0 0]; % Corrected to match rb_ACFmig
    rb_migTPT.TPT = -log(1-0.5) * [0 1 0 0]; % Corrected to match rb_TPT
    pb_migTPT.migrTPT = 0.8;
    Mb_migTPT = make_model(pb_migTPT, rb_migTPT, i, s, gps, prmb_migTPT.contmat);

    models = {M0, Mb, Mb_ACFdom, Mb_ACFmig, Mb_TPT, Mb_migTPT};
    prev_soln = init; % Initialize with initial conditions

    % Store M0’s end state for resetting
    M0_end_soln = []; % Temporary storage, reused each iteration

    for mi = 1:length(models)
        if mi == 1 % M0
            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, p0, [2024 2026], agg, prm0, sel, r0, false);
            [t, soln] = ode15s(geq, 2021:2041, prev_soln, opts);
            sdiff = diff(soln, [], 1);
            pops = sum(soln(:,1:i.nstates),2);
            incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
            prev_soln = soln(end, :);
            M0_end_soln = prev_soln; % Store M0’s end state
        else % Mb and subsequent models
            % Reset prev_soln to M0’s end state for Mb_ACFdom and beyond
            if mi >= 3
                prev_soln = M0_end_soln;
            end
            % Start from 2026
            if mi == 2
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb, [2026 2030], agg, prmb, sel, rb, true);
            elseif mi == 3
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb_ACFdom, [2026 2030], agg, prmb_ACFdom, sel, rb_ACFdom, true);
            elseif mi == 4
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb_ACFmig, [2026 2030], agg, prmb_ACFmig, sel, rb_ACFmig, true);
            elseif mi == 5
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb_TPT, [2026 2030], agg, prmb_TPT, sel, rb_TPT, true);
            elseif mi == 6
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb_migTPT, [2026 2030], agg, prmb_migTPT, sel, rb_migTPT, true);
            end
            [t, soln] = ode15s(geq, 2026:2041, prev_soln, opts);
            sdiff = diff(soln, [], 1);
            pops = sum(soln(:,1:i.nstates),2);
            inc2 = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1); % 15 values (2027–2041)
            % Pad with M0’s values for 2022–2026
            inc1 = incsto(1:5, ii, 1); % First 5 years from M0
            incsto(:, ii, mi) = [inc1; inc2];
            % Diagnostic check: ensure incidence not lower than Mb
            if mi >= 3 && any(incsto(6:end, ii, mi) < incsto(6:end, ii, 2))
                warning('Model %d incidence lower than Mb for parameter set %d', mi, ii);
            end
            prev_soln = soln(end, :);
        end
    end
end
fprintf('\n');

% Years vector based on incsto 
years = 2022:2041;

% Compute statistics
central_estimate = mean(incsto, 2);             
lowerbound = prctile(incsto, 2.5, 2);          
upperbound = prctile(incsto, 97.5, 2);        

% Plotting
ff = figure('Position', [577, 190, 1029, 732]); 
hold on;

colors = lines(2); % Only need 2 colors (M0 and Mb)
legendEntries = {'Baseline (M0)', 'New Action Plan (Mb)', 'Mb ACF Domestic', 'Mb ACF Migrant', 'Mb TPT Migrant', 'Mb Migrant Entry TPT'};

% Plot all models
for mi = 1:6
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    if mi <= 2 % M0 and Mb: solid lines with shaded areas
        lower = squeeze(lowerbound(:, 1, mi));
        upper = squeeze(upperbound(:, 1, mi));
        fill([years fliplr(years)], [lower' fliplr(upper')], ...
            colors(mi, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
        plot(years, central_estimate_model, 'LineWidth', 2, 'Color', colors(mi, :)); 
    else % Mb_ACFdom, Mb_ACFmig, Mb_TPT, Mb_migTPT: dashed lines, same color as Mb
        plot(years, central_estimate_model, 'LineWidth', 1.5, 'Color', colors(2, :), ...
             'LineStyle', '--'); 
    end
end

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
xlim([years(1) years(end)]);
ylim([0, 11]);
legend(legendEntries, 'FontWeight', 'bold', 'FontSize', 12, 'Location', 'best');

hold off;