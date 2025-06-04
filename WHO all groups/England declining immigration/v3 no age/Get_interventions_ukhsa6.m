clear all; load temp4.mat;

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
      
    init = aux.init;

    [p0,r0,prm0] = allocate_parameters(xx,p,r,xi,prm.scaling,prm);
    p0.prev_in_migr = 0;
    r0.gamma = r0.gamma_2015;
    r0.TPT = [0 r0.TPT2020rec 0];
    M0 = make_model2(p0, r0, i, s, gps, prm0.contmat);

    rb = r0; pb = p0; prmb = prm0;
    rb.ACF = -log(1-0.99) * [1 1 1 1];
    rb.TPT = -log(1-0.5) * [0 1 0 0];
    pb.migrTPT = 0.8;
    rb.muTx    = ( 109/4365) / ( 4044+103/4365) * rb.Tx;
    rb.ltfu    = ( 109/4365) / ( 4044+103/4365) * rb.Tx;   
    Mb = make_model2(pb, rb, i, s, gps, prmb.contmat);


            % Mb: treatment outcomes 
    rb_tx = r0; pb_tx = p0; prmb_tx = prm0;
    %rb_tx.TPT = -log(1-0.5) * [0 1 0 0]; 
    rb_tx.muTx    = ( 109/4365) / ( 4044+103/4365) * rb_tx.Tx;
    rb_tx.ltfu    = ( 109/4365) / ( 4044+103/4365) * rb_tx.Tx;   
    Mb_tx = make_model2(pb_tx, rb_tx, i, s, gps, prmb_tx.contmat);

        % Mb: and TPT migrant 
    rb_TPT = r0; pb_TPT = p0; prmb_TPT = prm0;
    rb_TPT.muTx    = ( 109/4365) / ( 4044+103/4365) * rb_TPT.Tx;
    rb_TPT.ltfu    = ( 109/4365) / ( 4044+103/4365) * rb_TPT.Tx;  
    rb_TPT.TPT = -log(1-0.5) * [0 1 0 0];
    Mb_TPT = make_model2(pb_TPT, rb_TPT, i, s, gps, prmb_TPT.contmat);


        % Mb: ACF migrant only
    rb_ACFmig = r0; pb_ACFmig = p0; prmb_ACFmig = prm0;
    rb_ACFmig.TPT = -log(1-0.5) * [0 1 0 0]; 
    rb_ACFmig.muTx    = ( 109/4365) / ( 4044+103/4365) * rb_ACFmig.Tx;
    rb_ACFmig.ltfu    = ( 109/4365) / ( 4044+103/4365) * rb_ACFmig.Tx;   
    rb_ACFmig.ACF = -log(1-0.99) * [0 1 1 1];
    Mb_ACFmig = make_model2(pb_ACFmig, rb_ACFmig, i, s, gps, prmb_ACFmig.contmat);

    % Mb: ACF domestic and migrant
    rb_ACFdom = r0; pb_ACFdom = p0; prmb_ACFdom = prm0;
    rb_ACFdom.TPT = -log(1-0.5) * [0 1 0 0]; 
    rb_ACFdom.muTx    = ( 109/4365) / ( 4044+103/4365) * rb_ACFdom.Tx;
    rb_ACFdom.ltfu    = ( 109/4365) / ( 4044+103/4365) * rb_ACFdom.Tx;   
    rb_ACFdom.ACF = -log(1-0.99) * [1 1 1 1];
    Mb_ACFdom = make_model2(pb_ACFdom, rb_ACFdom, i, s, gps, prmb_ACFdom.contmat);

    % Mb: and migrant entry TPT 
    rb_migTPT = r0; pb_migTPT = p0; prmb_migTPT = prm0;
    rb_migTPT.TPT = -log(1-0.5) * [0 1 0 0]; 
    rb_migTPT.muTx    = ( 109/4365) / ( 4044+103/4365) * rb_migTPT.Tx;
    rb_migTPT.ltfu    = ( 109/4365) / ( 4044+103/4365) * rb_migTPT.Tx;   
    rb_migTPT.ACF = -log(1-0.99) * [1 1 1 1];
    pb_migTPT.migrTPT = 0.8;
    Mb_migTPT = make_model2(pb_migTPT, rb_migTPT, i, s, gps, prmb_migTPT.contmat);


    models = {M0, Mb, Mb_tx, Mb_TPT, Mb_ACFmig, Mb_ACFdom};
    prev_soln = init;

    M0_end_soln = [];

    for mi = 1:length(models)
        if mi == 1 % M0
            geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, p0, [2024 2026], agg, prm0, sel, r0, false);
            [t, soln] = ode15s(geq, 2021:2041, prev_soln, opts);
            sdiff = diff(soln, [], 1);
            pops = sum(soln(:,1:i.nstates),2);
            incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
            mrtsto(:, ii, mi) = sdiff(:, i.aux.mort) * 1e5./pops(1:end-1);
            prev_soln = soln(end, :);
            M0_end_soln = prev_soln; 
        else 
            if mi >= 3
                prev_soln = M0_end_soln;
            end
            if mi == 2
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb, [2026 2030], agg, prmb, sel, rb, false);
            elseif mi == 3
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb_tx, [2026 2030], agg, prmb_tx , sel, rb_tx , false);
            elseif mi == 4
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb_TPT, [2026 2030], agg, prmb_TPT, sel, rb_TPT, false);
            elseif mi == 5
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb_ACFmig, [2026 2030], agg, prmb_ACFmig, sel, rb_ACFmig, false);
            elseif mi == 6
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb_ACFdom, [2026 2030], agg, prmb_ACFdom, sel, rb_ACFdom, false);
            elseif mi == 7
                geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, rin_vec, p0, pb_migTPT, [2026 2030], agg, prmb_migTPT, sel, rb_migTPT, false);
            end
            [t, soln] = ode15s(geq, 2026:2041, prev_soln, opts);
            sdiff = diff(soln, [], 1);
            pops = sum(soln(:,1:i.nstates),2);
            inc2 = sdiff(:, i.aux.inc(1)) * 1e5./pops(1:end-1);
            inc1 = incsto(1:5, ii, 1); %  M0 from 2022 to 2026
            mrt2 = sdiff(:, i.aux.mort) * 1e5./pops(1:end-1);
            mrt1 = mrtsto(1:5, ii, 1); %  M0 from 2022 to 2026
            incsto(:, ii, mi) = [inc1; inc2];
            mrtsto(:, ii, mi) = [mrt1; mrt2];
            prev_soln = soln(end, :);
        end
    end
end
fprintf('\n');


years = 2022:2041;
central_estimate = mean(incsto, 2);             
lowerbound = prctile(incsto, 2.5, 2);          
upperbound = prctile(incsto, 97.5, 2);        


ff = figure('Position', [577, 190, 1029, 732]); 
hold on;


colors = [
    0, 0, 1;    % M0
    1, 0, 0;    %  Mb
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
ylim([0, 11]);

hold off;

save intvn_res;