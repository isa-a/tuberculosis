% v2: Simulating using initial conditions starting from 2014, to ensure
% doing same simulations as in slidefig

clear all; load tmp0.mat; load Model_setup.mat;

% rin_vec  = [rin_vec, rin_vec(end)*5/9];
saving = 1;

obj = @(x) get_objective3(x, ref, prm, gps, prm.contmat, rin_vec, lhd);

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

midpt = false; 
if midpt
    xs = x3;
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
      
    init = aux.init;
    init = aux.soln(2010:2022==2014,:);

    [p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);
    p0.prev_in_migr = 0;
    r0.gamma = r0.gamma_2015;
    r0.TPT = [0 r0.TPT2020rec 0];
    M0 = make_model(p0, r0, i, s, gps, prm0.contmat);

    % Enhanced TPT, recent migrants
    r1 = r0; p1 = p0;
    r1.TPT = -log(1-0.25) * [0 1 0 0];
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
    r3.ACF = -log(1-0.25) * [0 1 1 1];
    M3 = make_model(p3, r3, i, s, gps, prm0.contmat);

    % ACF in domestic and foreign-born
    r4 = r3; p4 = p3;
    r4.ACF = -log(1-0.25) * [1 1 1 1];
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
        [t, soln] = ode15s(geq, 2014:2041, init, opts);

        sdiff = diff(soln, [], 1);
        pops = sum(soln(:,1:i.nstates),2);

        incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
        mrtsto(:, ii, mi) = sdiff(:, i.aux.mort) * 1e5./pops(1:end-1);
        prvsto(:, ii, mi) = sum(soln(:,s.allI),2) * 1e5./pops;

        if mi == 6
            incsource(ii,:) = sdiff(end,i.aux.incsources);
        end
    end

end
fprintf('\n');

% years = 2014:2040;
% prctile(incsto(years>=2014 & years<=2024,:,1),[2.5,50,97.5],2)
% return;
if saving 
    save intvn_resb;
end

years = 2015:2041;
central_estimate = mean(incsto, 2);             
lowerbound = prctile(incsto, 2.5, 2);          
upperbound = prctile(incsto, 97.5, 2);        

inds = find(years>=2014 & years<2024);

md = squeeze(central_estimate(inds, 1, 1))

vec1 = lowerbound(inds,1,1);
vec2 = upperbound(inds,1,1);

% [vec1, md, vec2]
% 
% return

% Extract prevalence estimates for Sharon, see email on 30 June 2025
prv_pct = prctile(prvsto,[2.5,50,97.5],2);
prv_tbl = [prv_pct(:,:,1), prv_pct(:,:,2), prv_pct(:,:,3), prv_pct(:,:,4), prv_pct(:,:,5), prv_pct(:,:,6)]';

incd_pct = prctile(incsto,[2.5,50,97.5],2);
incd_tbl = [incd_pct(:,:,1), incd_pct(:,:,2), incd_pct(:,:,3), incd_pct(:,:,4), incd_pct(:,:,5), incd_pct(:,:,6)]';

% Population numbers for Sharon
popn = [2022	57112542;
2023	57910211;
2024	58568219;
2025	59145411;
2026	59638662;
2027	60045068;
2028	60362009;
2029	60673451;
2030	60978896;
2031	61277976;
2032	61570538;
2033	61856818;
2034	62137307;
2035	62412168];

prevyrs = 2014:2021;
pop_old = [prevyrs; interp1(popn(:,1), popn(:,2), prevyrs, 'linear', 'extrap')]';

newyrs  = 2036:2040;
pop_new = [newyrs; interp1(popn(:,1), popn(:,2), newyrs, 'linear', 'extrap')]';

popn = [pop_old; popn; pop_new];
% figure; plot(popn(:,1), popn(:,2))

incd_nums = incd_tbl/1e5.*popn(:,2)'/1e3;
if saving
    writematrix(prv_tbl,   'prev_tbl.xlsx');
    writematrix(incd_tbl,  'incd_tbl.xlsx');
    writematrix(incd_nums, 'incd_nums.xlsx');
end


return;





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


