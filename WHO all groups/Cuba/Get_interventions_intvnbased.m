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
    nx  = 20;
    dx  = round(ix0/nx);
    xs  = xsto(ix0:dx:end,:);
end

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);

    init = aux.soln(end,:);

    [p0,r0] = allocate_parameters(xx,p,r,prm,xi);
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
%     rb = ra; pb = pa;
%     pb.migrTPT = 0.25;
%     Mb = make_model(pb, rb, i, s, gps, prm.contmat);

    % TPT in vulnerables
%     rc = rb; pc = pb;
%     rc.TPT = TPTcov * [0 1 0 1];
%     Mc = make_model(pc, rc, i, s, gps, prm.contmat);

    %diag delay
    re = ra; pe = pa; 
    re.ACF = -log(1-0.99) * [1 1 1 1];
%     re.TPT = TPTcov * [0 1 0 1];
%     pe.migrTPT = 0.25; 
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
%     rg.ACF = -log(1-0.99) * [1 1 1 1];
%     rg.TPT = -log(1-0.5) * [1 1 1 1];
%     pg.migrTPT = 0.8;
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
legendEntries = {'Baseline', 'Oman NSP', 'Expanded deployment of current tools', 'New TPT and treatments'};
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
    %ylim([0, 10]);
    yline(1,'k--','LineWidth', 2, 'HandleVisibility', 'off');
    yline(0.1,'k--','LineWidth', 2, 'HandleVisibility', 'off');
    
    pause;  % Waits for the user to press a key
end

hold off;



return;


incmat   = permute(prctile(incsto,[2.5,50,97.5],2),[2,1,3]);
incmat2  = permute(prctile(incsto2,[2.5,50,97.5],2),[2,1,3]);
incmatRR = permute(prctile(incstoRR,[2.5,50,97.5],2),[2,1,3]);

mrtmat = permute(prctile(mrtsto,[2.5,50,97.5],2),[2,1,3]);


% -------------------------------------------------------------------------
% --- Plot figure of incidence and mortality impacts ----------------------

ff=figure('Position', [577,   190 ,   1029 ,732]); lw = 1.5; fs = 14;
% allmat = cat(4,incmat,incmat2);
allmat = cat(4,incmat,incmat2,incmatRR);

cols = linspecer(size(allmat,3));
cols1 = linspecer(size(mrtmat, 3));
xx = [2022:2040];

lgs = {'Baseline (continuation of existing TB services)','TPT, in-country migrants, 90% annually','+ TPT in 100% of migrants, pre-arrival','+ TPT in contacts','+ TPT in vulnerable population, 90% annually','+ Find and treat active TB in migrants and vulnerable', '+ TPT in general population, 90% annually'};


plotopts = {'All incidence','RR incidence','Alternative incidence'};
plotopt  = plotopts{1};




if strcmp(plotopt, 'All incidence')
    
    % --- Single incidence plot, showing one by one
    hold on;
    for ii = 1:size(allmat,3)
        plt = allmat(:,:,ii,1);
        lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
        %jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
        
        ylim([0 7.5]);
        xlim([2022 2040]);
        set(gca,'fontsize',fs);
        
        yline(1,'k--'); yline(0.1,'k--');
        ylabel('Incidence rate per 100,000 population');
        xlabel('Year');
        
        legend(lg,lgs(1:ii),'Location','Southwest');
        pause;
    end
    
elseif strcmp(plotopt, 'RR incidence')
    
    % --- Incidence of RR-TB
    hold on;
    
    selinds = [1:7];
    allplt = squeeze(allmat(:,:,selinds,3));
    lgsel  = lgs(selinds);

    
    cols(1,:) = [0.545, 0, 0]; % Dark Red for Baseline
    cols(2,:) = [0, 0, 1]; % Blue for TPT, in-country migrants, 90% annually
    cols(3,:) = [0, 0.75, 0]; % GREEN
    
    for ii = 1:size(allplt,3)
        plt = allplt(:,:,ii);
        if ii == 2
            lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
        else
            lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
        end
        
        xlim([2022 2040]);
        set(gca,'fontsize',fs);
        ylim([0.08 0.125])
        yline(0.1,'k--');
        ylabel('Incidence rate per 100,000 population');
        %title('Additional coverage of TPT', 'FontWeight', 'bold');
        xlabel('Year', 'FontWeight', 'bold');
        set(gca, 'FontWeight', 'bold');
        
                
        legend(lg,lgsel(1:ii),'Location','NorthEast');
    end
    
elseif strcmp(plotopt, 'Alternative incidence')
    
    % --- Two ways of counting
    tis = {'All incidence','Incidence without treatment history'};
    for is = 1:2
        subplot(1,2,is); hold on;
    
        for ii = 1:size(allmat,3)
            plt = allmat(:,:,ii,is);
            lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
            jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
        end
        yl = ylim; yl(1) = 0; ylim(yl); ylims(is,:) = ylim;
        xlim([2022 2040]);
        set(gca,'fontsize',fs);
    
        title(tis{is});
        
        yline(1,'k--');
        yline(0.1,'k--');
    
    end
    subplot(1,2,2);
    ylim(ylims(1,:));
    
    % legend(lg, 'Baseline','TPT, recent migrants','+ Case-finding, active TB','+ TPT, new migrants (hypothetical)','+ TPT, domestic (hypothetical)', 'Elimination target','location','SouthWest');
    legend(lg,'Baseline','TPT, in-country migrants, 50% annually','+ TPT in 100% of migrants, pre-arrival','+ TPT in vulnerable population, 50% annually','+ ACF in migrants and vulnerable', '+ TPT in general population, 50% annually');
    subplot(1,2,1);
    ylabel('Rate per 100,000 population');

end













return;

tis = {'All incidence','Incidence without treatment history'};
for is = 1:2
    subplot(1,2,is); hold on;

    for ii = 1:size(allmat,3)
        plt = allmat(:,:,ii,is);
        lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
        jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
    end
    yl = ylim; yl(1) = 0; ylim(yl); ylims(is,:) = ylim;
    xlim([2022 2040]);
    set(gca,'fontsize',fs);

    title(tis{is});

    yline(1,'k--');
    yline(0.1,'k--');

end
subplot(1,2,2);
ylim(ylims(1,:));

% legend(lg, 'Baseline','TPT, recent migrants','+ Case-finding, active TB','+ TPT, new migrants (hypothetical)','+ TPT, domestic (hypothetical)', 'Elimination target','location','SouthWest');
legend(lg,'Baseline','TPT, in-country migrants, 50% annually','+ TPT in 100% of migrants, pre-arrival','+ TPT in vulnerable population, 50% annually','+ ACF in migrants and vulnerable', '+ TPT in general population, 50% annually');
subplot(1,2,1);
ylabel('Rate per 100,000 population');


% --- Find remaining sources of incidence
vec = squeeze(props(end,:,end));
tmp  = abs(vec)/sum(abs(vec));
tmp2 = sortrows([tmp; 1:length(tmp)]',-1);

lbls = {}; ind = 1;
set1 = {'ds','rr'};
set2 = {'dom','migr','vuln'};
set3 = {'L','P','R','T'};
for is3 = 1:length(set3)
    for is2 = 1:length(set2)
        for is1 = 1:length(set1)
            lbls{ind} = [set1{is1}, ' ', set2{is2}, ' ', set3{is3}];
            ind = ind+1;
        end
    end
end

inds = tmp2(1:6,2)'
tmp2(1:6,1)'
lbls(inds)





return;

% -------------------------------------------------------------------------
% --- Plot figure of incidence components as of 2030 ----------------------

tmp1 = reshape(props(:,:,end),3,5);                                        % Dims: 1.Dom/migr_rect/migr_long, 2.Lf,Pf,Ls,Ps,R
tmp2 = [tmp1(1,:); sum([tmp1(2,:); tmp1(3,:)],1)];                         % Dims: 1.Dom/all migr, 2.Lf,Pf,Ls,Ps,R
tmp3 = [sum(tmp2(:,[1,3]),2), sum(tmp2(:,[2,4]),2), tmp2(:,end)];          % Dims: 1.Dom/all migr, 2.All L, All P, R
tmp4 = [tmp3(1,:), tmp3(2,:)];
labels = {'UK-born without treatment history', 'UK-born after TPT', 'UK-born after TB Rx', 'Migrants without treatment history', 'Migrants after TPT', 'Migrants after TB Rx'};
figure; pie(tmp4);
legend(labels,'Location','NorthWest','Orientation','vertical');
title('Sources of incidence in 2035 with all interventions combined')