clear all; load optim_res1000_2.mat;

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
% Initialize storage for additional interventions
incsto = zeros(20, size(xs,1), 6); % 3 original models + 3 additional interventions (1 for Ma, 2 for Mb)
for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);
      
    init  = aux.init;
    incd  = aux.incd;
    incdt = aux.incdt;

    [p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);
    p0.prev_in_migr = 0;
    r0.gamma        = r0.gamma_2015;
    r0.TPT          = [0 r0.TPT2020rec 0];
    M0              = make_model(p0, r0, i, s, gps, prm0.contmat);

    TPTcov = -log(1-0.25);
    ACFcov = -log(1-0.50); 

    % UKHSA action plan (from now)
    ra = r0; pa = p0; prma = prm0;
    ra.TPT = TPTcov * [0 1 0 0];
    ra.ACF = ACFcov * [1 0 0 0];
    Ma = make_model(pa, ra, i, s, gps, prma.contmat);

    % Expanded deployment of current tools (2027 – 2030)
    rb = ra; pb = pa; prmb = prma;
    rb.ACF = -log(1-0.99) * [1 1 1 1];
    rb.TPT = -log(1-0.5) * [1 1 1 1];
    pb.migrTPT = 0.8;
    Mb = make_model(pb, rb, i, s, gps, prmb.contmat);

    % Define additional models for individual interventions
    % For Ma: TPT only
    ra_TPT = r0; pa_TPT = p0; prma_TPT = prm0;
    ra_TPT.TPT = TPTcov * [0 1 0 0];
    Ma_TPT = make_model(pa_TPT, ra_TPT, i, s, gps, prma_TPT.contmat);

    % For Mb: ACF only, then TPT only
    rb_ACF = ra; pb_ACF = pa; prmb_ACF = prma;
    rb_ACF.ACF = -log(1-0.99) * [1 1 1 1];
    Mb_ACF = make_model(pb_ACF, rb_ACF, i, s, gps, prmb_ACF.contmat);

    rb_TPT = rb_ACF; pb_TPT = pb_ACF; prmb_TPT = prmb_ACF;
    rb_TPT.TPT = -log(1-0.5) * [1 1 1 1];
    Mb_TPT = make_model(pb_TPT, rb_TPT, i, s, gps, prmb_TPT.contmat);

    models = {M0, Ma, Mb, Ma_TPT, Mb_ACF, Mb_TPT};

    for mi = 1:length(models)
        if mi <= 3 % Original models
            if mi < 3
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pa, [2024 2027], agg, prm0, sel, r0, false);
                [t, soln] = ode15s(geq, 2021:2041, init, opts);
                sdiff = diff(soln, [], 1);
                pops = sum(soln(:,1:i.nstates),2);
                incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
            else
                % Expanded deployment (mi=3): run Ma to 2027 then Mb to 2041
                geq1 = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{2}, rin_vec, p0, pa, [2024 2027], agg, prma, sel, ra, false);
                [~, soln1] = ode15s(geq1, 2021:2027, init, opts);
                sdiff_ma = diff(soln1, [], 1);
                pops1 = sum(soln1(:, 1:i.nstates), 2);
                inc1 = sdiff_ma(:, i.aux.inc(1)) * 1e5 ./ pops1(1:end-1);
        
                geq2 = @(t,in) goveqs_scaleup(t, in, i, s, Ma, models{3}, rin_vec, pa, pb, [2027 2030], agg, prmb, sel, rb, true);
                [~, soln2] = ode15s(geq2, 2027:2041, soln1(end, :), opts);
                sdiff_mb = diff(soln2(:, i.aux.inc(1)), [], 1);
                pops2 = sum(soln2(:, 1:i.nstates), 2);
                inc2 = sdiff_mb * 1e5 ./ pops2(1:end-1);
        
                incsto(:, ii, 3) = [inc1; inc2];
            end
        elseif mi == 4 % Ma_TPT
            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pa_TPT, [2024 2027], agg, prm0, sel, r0, false);
            [t, soln] = ode15s(geq, 2021:2041, init, opts);
            sdiff = diff(soln, [], 1);
            pops = sum(soln(:,1:i.nstates),2);
            incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
        elseif mi >= 5 % Mb_ACF and Mb_TPT
            % Run Ma to 2027 first
            geq1 = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{2}, rin_vec, p0, pa, [2024 2027], agg, prma, sel, ra, false);
            [t1, soln1] = ode15s(geq1, 2021:2027, init, opts);
            sdiff_ma = diff(soln1, [], 1);
            pops1 = sum(soln1(:, 1:i.nstates), 2);
            inc1 = sdiff_ma(:, i.aux.inc(1)) * 1e5 ./ pops1(1:end-1); % 6 values (2022–2027)
            
            % Run Mb_ACF or Mb_TPT from 2027 to 2041
            if mi == 5 % Mb_ACF
                geq2 = @(t,in) goveqs_scaleup(t, in, i, s, Ma, models{mi}, rin_vec, pa, pb_ACF, [2027 2030], agg, prmb_ACF, sel, rb_ACF, true);
                [t2, soln2] = ode15s(geq2, 2027:2041, soln1(end, :), opts);
                sdiff_mb = diff(soln2(:, i.aux.inc(1)), [], 1);
                pops2 = sum(soln2(:, 1:i.nstates), 2);
                inc2 = sdiff_mb * 1e5 ./ pops2(1:end-1); % 14 values (2028–2041)
                incsto(:, ii, mi) = [inc1; inc2]; % 6 + 14 = 20 values
            else % mi == 6, Mb_TPT
                % Run Mb_ACF to 2030 to get the starting point
                geq_acf = @(t,in) goveqs_scaleup(t, in, i, s, Ma, models{5}, rin_vec, pa, pb_ACF, [2027 2030], agg, prmb_ACF, sel, rb_ACF, true);
                [~, soln_acf] = ode15s(geq_acf, 2027:2030, soln1(end, :), opts);
                % Run Mb_TPT from 2027 to 2041, starting from Ma's end
                geq2 = @(t,in) goveqs_scaleup(t, in, i, s, Ma, models{6}, rin_vec, pa, pb_TPT, [2027 2030], agg, prmb_TPT, sel, rb_TPT, true);
                [t2, soln2] = ode15s(geq2, 2027:2041, soln1(end, :), opts);
                sdiff_mb = diff(soln2(:, i.aux.inc(1)), [], 1);
                pops2 = sum(soln2(:, 1:i.nstates), 2);
                inc2 = sdiff_mb * 1e5 ./ pops2(1:end-1); % 14 values (2028–2041)
                incsto(:, ii, mi) = [inc1; inc2]; % 6 + 14 = 20 values
            end
        end
    end
end
fprintf('\n');

% Years vector based on incsto 
years = 2022:(2022 + size(incsto, 1) - 1);

% Compute statistics
central_estimate = mean(incsto, 2);             
lowerbound = prctile(incsto, 2.5, 2);          
upperbound = prctile(incsto, 97.5, 2);        

% Plotting
ff = figure('Position', [577, 190, 1029, 732]); 
hold on;

colors = lines(3);
legendEntries = {'Baseline', 'UKHSA Action Plan', 'New Action Plan', ...
                 'Ma TPT only', 'Mb ACF only', 'Mb TPT only'};

for mi = 1:6
    % Extract time series for current model
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    if mi <= 3
        % Original models: plot shaded areas and solid lines
        lower = squeeze(lowerbound(:, 1, mi));
        upper = squeeze(upperbound(:, 1, mi));
        fill([years fliplr(years)], [lower' fliplr(upper')], ...
            colors(mi, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
        plot(years, central_estimate_model, 'LineWidth', 2, 'Color', colors(mi, :)); 
    else
        % Intervention models: plot dashed lines only
        % Use color of parent model (Ma for mi=4, Mb for mi=5,6)
        color_idx = 2 + (mi >= 5); % Ma (2) for mi=4, Mb (3) for mi=5,6
        plot(years, central_estimate_model, 'LineWidth', 1.5, 'Color', colors(color_idx, :), ...
             'LineStyle', '--'); 
    end
end

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
xlim([years(1) years(end)]);
ylim([0, 11]);
legend(legendEntries, 'FontWeight', 'bold', 'FontSize', 12, 'Location', 'best');

hold off;