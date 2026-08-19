clear all; load tmp0.mat;

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
incsto = nan(20, size(xs,1), 6);
mrtsto = nan(20, size(xs,1), 6);
mk = round(size(xs,1)/25);
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

    % Enhanced TPT, recent migrants
    r1 = r0; p1 = p0;
    r1.TPT = -log(1-0.5) * [0 1 0 0];
    M1 = make_model(p1, r1, i, s, gps, prm0.contmat);

    % Improved Tx outcomes
    r2 = r1; p2 = p1;

    % Find new treatment outcomes
    vec = [r1.Tx, r1.ltfu, r1.muTx]; props = vec/sum(vec);

    props(end) = 0.01; props(2) = props(2)/2; props(1) = 1 - sum(props(2:3));
    newrates = r1.Tx*props/props(1);
    r2.ltfu = newrates(2);
    r2.muTx = newrates(3);
    M2 = make_model(p2, r2, i, s, gps, prm0.contmat);

    % ACF in foreign-born
    r3 = r2; p3 = p2;
    r3.ACF = -log(1-0.99) * [0 1 1 1];
    M3 = make_model(p3, r3, i, s, gps, prm0.contmat);

    % ACF in domestic and foreign-born
    r4 = r3; p4 = p3;
    r4.ACF = -log(1-0.99) * [1 1 1 1];
    M4 = make_model(p4, r4, i, s, gps, prm0.contmat);

    % Pre-entry TPT
    r5 = r4; p5 = p4;
    p5.migrTPT = 0.6;
    M5 = make_model(p5, r5, i, s, gps, prm0.contmat);

    models = {M0, M1, M2, M3, M4, M5};
    % prev_soln = init;

    M0_end_soln = [];

    for mi = 1:length(models)

        geq = @(t,in) goveqs_scaleupb(t, in, i, s, M0, models{mi}, rin_vec, [2026 2030], agg, prm0, sel, r0, p0, false);
        [t, soln] = ode15s(geq, 2021:2041, init, opts);

        sdiff = diff(soln, [], 1);
        pops = sum(soln(:,1:i.nstates),2);

        incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
        mrtsto(:, ii, mi) = sdiff(:, i.aux.mort) * 1e5./pops(1:end-1);
    end

end
fprintf('\n');
save intvn_resb;

years = 2022:2041;
central_estimate = mean(incsto, 2);             
lowerbound = prctile(incsto, 2.5, 2);          
upperbound = prctile(incsto, 97.5, 2);        


ff = figure('Position', [577, 190, 1029, 732]); 
hold on;


colors = [
    0, 0, 1;    % M0
    0.3, 0, 0;    %  Mb
    0, 0.4470, 0.7410;    %  Mb_TPT
    0, 0, 0;    %  Mb_tx
    1, 0, 1;    %  Mb_ACFmig
    0.9290 0.6940 0.1250;    %  Mb_ACFdom
    0.5, 0, 0.5 %  Mb_migTPT
];



for mi = 1:6
    central_estimate_model = squeeze(central_estimate(:, 1, mi));
    if mi <= 6
        lower = squeeze(lowerbound(:, 1, mi));
        upper = squeeze(upperbound(:, 1, mi));
        fill([years fliplr(years)], [lower' fliplr(upper')], ...
            colors(mi, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
        plot(years, central_estimate_model, 'LineWidth', 2, 'Color', colors(mi, :)); 
    else 
        plot(years, central_estimate_model, 'LineWidth', 1.5, 'Color', colors(mi, :), ...
             'LineStyle', '-'); 
    end
end

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
xlim([years(1) years(end)]);
% ylim([0, 11]);

% hold off;

years = [2015, 2016,2017, 2018,2019, 2020, 2021, 2022, 2023, 2024];
incs1 = [5734,5621,5066,4610,4704,4124,4407,4375,4855,0] * 0.001538;
incs1(end)   = 9.5; 

plot(years, incs1, '.', 'MarkerSize', 20);


