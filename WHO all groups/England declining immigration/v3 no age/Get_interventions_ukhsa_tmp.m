load optim_res1000.mat;

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
      
    init    = aux.soln(2022 - 2010 + 1, :);

    [p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);
    r0.gamma = r0.gamma_2015;
    p0.LTBIdec = 0;
    p0.prev_in_migr = 0.003;
    M0 = make_model(p0,r0,i,s,gps,prm0.contmat);

    % TPTcov = -log(1-0.25);
    % ACFcov = -log(1-0.50); 
    % 
    % % UKHSA action plan (from now)
    % % 25% annually of recent migrants (arrived within the last 5 years)
    % % Reduce the average delay to diagnosis from 75 days to 56 days
    % ra = r0; pa = p0; prma = prm0;
    % ra.TPT = TPTcov * [0 1 0 0];
    % ra.ACF = ACFcov * [1 0 0 0];
    % Ma = make_model(pa, ra, i, s, gps, prma.contmat);
    % 
    % %Expanded deployment of current tools (2027 – 2030)
    % %Scale up diagnosis efforts to detect 100% of cases in all population groups
    % %50% annually in all TB infected: UK born, Non-UK born
    % %Pre-entry screening to cover 80% of new migrants (80% with LTBI completing treatment)
    % rb = ra; pb = pa; prmb = prma;
    % rb.ACF = -log(1-0.99) * [1 1 1 1];
    % rb.TPT = -log(1-0.5) * [1 1 1 1];
    % pb.migrTPT = 0.8;
    % Mb = make_model(pb, rb, i, s, gps, prmb.contmat);

    models = {M0};

    for mi = 1:length(models)
        if mi < 3

            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pa, [2024 2027], agg, prm0, sel, r0, false);
            [t, soln] = ode15s(geq, 2022:2041, init, opts);
            sdiff = diff(soln, [], 1);
            incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5;
        else
            % --- expanded deployment (mi=3): run Ma to 2027 then Mb to 2041
            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{2}, rin_vec, p0, pa, [2024 2027], agg, prma, sel, ra, false);
            [~, soln] = ode15s(geq, 2022:2027, init, opts);
            sdiff_ma = diff(soln, [], 1);
            
            geq = @(t,in) goveqs_scaleup(t,in,i,s,M0,models{3},rin_vec,p0,pb,[2027 2030],agg,prmb,sel,rb,true);
            [~, soln2] = ode15s(geq, 2027:2041, soln(end,:), opts);
            sdiff_mb = diff(soln2(:, i.aux.inc(1)));

            incsto(:, ii, 3) = [sdiff_ma(:, i.aux.inc(1)); sdiff_mb] * 1e5;
        end
    end
end
fprintf('\n');


mat  = prctile(incsto,[2.5,50,97.5],2);
mat2 = squeeze(mat(1,2,1))


return;


%  years vector based on incsto 
years = 2022:(2022 + size(incsto, 1) - 1);


central_estimate = mean(incsto, 2);             
lowerbound = prctile(incsto, 2.5, 2);          
upperbound = prctile(incsto, 97.5, 2);        


ff = figure('Position', [577, 190, 1029, 732]); 
hold on;


colors = lines(3);
legendEntries = {'Baseline', 'UKHSA Action Plan', 'New Action Plan'};


for mi = 1:3
    % extract time series for current model
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    lower = squeeze(lowerbound(:, 1, mi));
    upper = squeeze(upperbound(:, 1, mi));

    %  shaded area for percentiles
    fill([years fliplr(years)], [lower' fliplr(upper')], ...
        colors(mi, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
    
    %  central estimate
    plot(years, central_estimate_model, 'LineWidth', 2, 'Color', colors(mi, :)); 
end


xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
xlim([years(1) years(end)]);
ylim([0, 11]);
legend(legendEntries, 'FontWeight', 'bold', 'FontSize', 12);

hold off;