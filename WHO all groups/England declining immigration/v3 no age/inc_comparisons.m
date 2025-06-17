obj = @(x) get_objective3(x, ref, prm, gps, prm.contmat, rin_vec, lhd);

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

midpt = false; 
if midpt
    xs = x3;
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
      
    init  = aux.init;
    incd  = aux.incd;
    incdt = aux.incdt;

    [p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);
    p0.prev_in_migr = 0;
    r0.gamma        = r0.gamma_2015;
    r0.TPT          = [0 r0.TPT2020rec 0];
    M0              = make_model(p0, r0, i, s, gps, prm0.contmat);

    models = {M0};

    for mi = 1:length(models)
        if mi == 1 
            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, p0, [2024 2027], agg, prm0, sel, r0, false);
            [t, soln] = ode15s(geq, 2014:2024, init, opts);
            sdiff = diff(soln, [], 1);
            pops = sum(soln(:,1:i.nstates),2);
            incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
        end
    end
end
fprintf('\n');


all_years = 2015:(2015 + size(incsto, 1) - 1);


central_estimate = mean(incsto, 2);             
lowerbound = prctile(incsto, 2.5, 2);           
upperbound = prctile(incsto, 97.5, 2);          


full_incs1 = [5734,5621,5066,4610,4704,4124,4407,4375,4855,0] * 0.001538;
full_incs1(end)   = 9.5; 


start_year = 2015;
idx = find(all_years >= start_year);
years = all_years(idx);

central_estimate = central_estimate(idx,:,:);
lowerbound = lowerbound(idx,:,:);
upperbound = upperbound(idx,:,:);
incs1 = full_incs1(idx);


ff = figure('Position', [577, 190, 1029, 732]); 
hold on;

colors = lines(4);

for mi = 1:1
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    if mi <= 3
        lower = squeeze(lowerbound(:, 1, mi));
        upper = squeeze(upperbound(:, 1, mi));
        fill([years fliplr(years)], [lower' fliplr(upper')], ...
            colors(mi, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
        plot(years, central_estimate_model, 'LineWidth', 2, 'Color', colors(mi, :)); 
    else
        color_idx = 2 + (mi >= 5);
        plot(years, central_estimate_model, 'LineWidth', 1.5, 'Color', colors(color_idx, :), ...
             'LineStyle', '--'); 
    end
end

%  data
plot(years, incs1, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');

% previois inc
x3data = load('x3_output.mat');  
x3_t = x3data.central_estimate_model(x3data.years >= start_year);
x3_y = x3data.years(x3data.years >= start_year);


plot(x3_y, x3_t, '-', ...
     'LineWidth', 2, 'Color', [1 0 0], ...
     'DisplayName', 'x3 from saved file');


xlabel('Year', 'FontWeight', 'bold', 'FontSize', 18);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 18);
xlim([years(1) years(end)]);
ylim([0, 10]);

legend({'Incidence modelled using immigration based on visas issued', 'Notification-based incidence', 'Incidence modelled using House of Commons migration data'}, ...
       'Location', 'NorthEast', 'FontSize', 18);
set(gca, 'FontSize', 12);  
hold off;

return;
