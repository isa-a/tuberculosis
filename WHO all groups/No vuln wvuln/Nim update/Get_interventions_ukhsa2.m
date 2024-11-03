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

    %diag delay
    re = rc; pe = pc; 
    re.ACF = ACFcov * [1 0 0 0];
    re.TPT = TPTcov * [0 1 0 1];
    pe.migrTPT = 0.25; 
    Me = make_model(pe, re, i, s, gps, prm.contmat);

    % 2027
    rf = re; pf = pe;
    rf.ACF = -log(1-0.99) * [1 1 1 1];
    rf.TPT = -log(1-0.5) * [1 1 1 1];
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


    models = {M0, Me, Mf, Mg};    
    
    for mi = 1:length(models)
        if mi == length(models)
            % Run Mf to 2030, then Mg
            geq_mf = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi-1}, p0, pf, [2024 2030], agg, sel, r0);
            [t1, soln1] = ode15s(geq_mf, 2022:2030, init, opts);

            % Final solution for Mf
            init_final = soln1(end,:);

            % Mg from 2030 
            geq_mg = @(t,in) goveqs_scaleup(t, in, i, s, models{mi-1}, models{mi}, p0, pg, [2030 2033], agg, sel, r0);
            [t2, soln2] = ode15s(geq_mg, 2030:2041, init_final, opts);

            % Combine
            t = [t1; t2(2:end)];
            soln = [soln1; soln2(2:end,:)];
        elseif mi == length(models)-1
            % For Mf, run Me up to 2027, then Mf 
            geq_me = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi-1}, p0, pe, [2024 2027], agg, sel, r0);
            [t1, soln1] = ode15s(geq_me, 2022:2027, init, opts);

            % Final solution for Me
            init_final = soln1(end,:);

            % Mf from 2027
            geq_mf = @(t,in) goveqs_scaleup(t, in, i, s, models{mi-1}, models{mi}, p0, pf, [2027 2030], agg, sel, r0);
            [t2, soln2] = ode15s(geq_mf, 2027:2041, init_final, opts);

            % Combine
            t = [t1; t2(2:end,:)];
            soln = [soln1; soln2(2:end,:)];
        else
            % All others
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


years = 2022:(2022 + size(incsto, 1) - 1);

% Percentiles
central_estimate = mean(incsto, 2);              
lower_bound = prctile(incsto, 2.5, 2);          
upper_bound = prctile(incsto, 97.5, 2);          

% Plot
ff = figure('Position', [577, 190, 1029, 732]); 
hold on;


colors = lines(length(models));
legendEntries = {'Baseline', 'UKHSA Action Plan', 'Max scaleup of new tools', 'New TPT and treatments'};
plotHandles = gobjects(length(models),1);


for mi = 1:length(models)
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    lower_bound_model = squeeze(lower_bound(:, 1, mi));
    upper_bound_model = squeeze(upper_bound(:, 1, mi));
    
    if mi == length(models)
        idx_2030 = find(years == 2030);
        
        central_estimate_model(1:idx_2030) = central_estimate(1:idx_2030, 1, mi-1);
        lower_bound_model(1:idx_2030) = lower_bound(1:idx_2030, 1, mi-1);
        upper_bound_model(1:idx_2030) = upper_bound(1:idx_2030, 1, mi-1);
    end

    fillHandle = fill([years fliplr(years)], [lower_bound_model' fliplr(upper_bound_model')], ...
        colors(mi, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
    

    plotHandle = plot(years, central_estimate_model, 'LineWidth', 2, 'Color', colors(mi, :)); 
    

    plotHandles(mi) = plotHandle;
    legend(plotHandles(1:mi), legendEntries(1:mi), 'FontWeight', 'bold', 'FontSize', 12);
    

    xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
    xlim([years(1) years(end)]);
    ylim([0, 9]);
    yline(1,'k--','LineWidth', 2, 'HandleVisibility', 'off');
    yline(0.1,'k--','LineWidth', 2, 'HandleVisibility', 'off');
    
    pause;  % Waits for the user to press a key
end

hold off;