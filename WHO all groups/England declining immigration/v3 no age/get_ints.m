clear all; load optim_res_MAIN.mat; load Model_setup;

obj = @(x) get_objective2(x, ref, prm, gps, prm.contmat, lhd);

% set1 = {'ds','rr'};
set1 = {'ds'};
set2 = {'dom','migr','vuln'};
set3 = {'L','P','R','T'};
[inci, incs, incd, lim] = get_addresses({set3, set2, set1}, [], [], [], 0);
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

midpt = true; 
if midpt
    % inds = find(outsto==max(outsto));
    % xs = xsto(inds(1),:);
    xs = x0sto(2,:);
else
    ix0 = size(xsto,1)/2;
    nx  = 200;
    dx  = round(ix0/nx);
    xs  = xsto(ix0:dx:end,:);
end

for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);
    

    init = aux.soln(end,:);

    [p0,r0] = allocate_parameters(xx,p,r,xi);
    r0.gamma = r0.gamma_2020;
    M0 = make_model(p0,r0,i,s,gps,prm.contmat);
    
    % ---------------------------------------------------------------------
    % --- Model intervention
    
    TPTcov = -log(1-0.99); %TPTcov = 100;
    ACFcov = -log(1-0.99); %ACFcov = 100;

    % social protection
    ra1 = r0; pa1 = p0;
    ra1.TPT = TPTcov*[0 0 0 1];
    pa1.TPTeff = 0.8;
    Ma1 = make_model(pa1,ra1,i,s,gps,prm.contmat);    

    % TPT in recent migrants
    ra = ra1; pa = pa1;
    ra.TPT = TPTcov*[0 1 0 0];
    Ma = make_model(pa,ra,i,s,gps,prm.contmat);
    
    % TPT at point of entry
    rb = ra; pb = pa;
    pb.migrTPT = 1;
    Mb = make_model(pb,rb,i,s,gps,prm.contmat);

    % TPT in longer-term migrants
    rb1 = rb; pb1 = pb;
    pb1.TPT = TPTcov*[0 1 1 1];
    Mb1 = make_model(pb1,rb1,i,s,gps,prm.contmat);

    % TPT in contacts
    rb2 = rb1; pb2 = pb1;
    rb2.progression  = rb2.progression*0.9;
    rb2.reactivation = rb2.reactivation*0.9;
    Mb2 = make_model(pb2,rb2,i,s,gps,prm.contmat);

    % TPT in vulnerables
    rc = rb2; pc = pb2;
    rc.TPT = TPTcov*[0 1 1 1];
    Mc = make_model(pc,rc,i,s,gps,prm.contmat);

    % ACF in migrants and vulnerables
    rd = rc; pd = pc;
    rd.ACF  = ACFcov*[0 1 0 1];
    rd.ACF2 = rd.ACF;
    Md = make_model(pd,rd,i,s,gps,prm.contmat);

    % ACF in general population
    re = rd; pe = pd;
    re.ACF  = ACFcov*[1 1 1 1];
    re.ACF2 = rd.ACF;
    re.TPT  = TPTcov*[1 1 1 1];
    Me = make_model(pe,re,i,s,gps,prm.contmat);

    % new tools, lower reactivation progression rates, higher TPT eff
    % vaccine - f = f*(1-cov*eff)
    % reduced relapse
    rf = re; pf = pe;
    pf.TPTeff = 0.8;
%     rf.relapse = re.relapse/2/2/2;
    rf.progression  = rf.progression*(1 - (0.9*0.9));
    rf.reactivation = rf.reactivation*(1 - (0.9*0.9));
    Mf = make_model(pf,rf,i,s,gps,prm.contmat);
    

    models = {M0, Ma1, Ma, Mb, Mb2, Mc, Md, Me, Mf};    
    for mi = 1:length(models)
    if mi == length(models) 
        % run Me
        geq_me = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi-1}, p0, pa, [2024 2029], agg, sel, r0);
        [t1, soln1] = ode15s(geq_me, 2022:2028, init, opts);

        % final solution for Me
        init_final = soln1(end,:);

        % mf 2027 onwards
        geq_mf = @(t,in) goveqs_scaleup(t, in, i, s, models{mi-1}, models{mi}, p0, pa, [2028 2032], agg, sel, r0);
        [t2, soln2] = ode15s(geq_mf, 2028:2041, init_final, opts);

        t = [t1; t2(2:end)];    % combine time from 2022 to 2028 with 2029 onwards
        soln = [soln1; soln2(2:end,:)];     %combine the two solutions
    else
        % every other scenario
        geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, pa, [2024 2029], agg, sel, r0);
        [t, soln] = ode15s(geq, 2022:2041, init, opts);
    end

    endsolsto(mi,:) = soln(end,:);

    sdiff = diff(soln, [], 1);
    incsto(:, ii, mi) = sdiff(:, i.aux.inc(1)) * 1e5;
    
    mrtsto(:, ii, mi) = sdiff(:, i.aux.mort) * 1e5;
    
    % Get proportions from different sources
    vec = sdiff(end, i.aux.incsources) * 1e5;
    props(ii, :, mi) = vec / sum(vec);
    end
end
fprintf('\n');
return;
ff=figure('Position', [577,   190 ,   1029 ,732]); lw = 1.5; fs = 14;
%figure; 
plot(squeeze(incsto), 'LineWidth', 2); % Make curves thicker/bolder
years = 2022:2035;
yl = ylim; 
yl(1) = 0; 
ylim(yl);

% set(gca, 'XTick', 1:size(squeeze(incsto), 1));
% set(gca, 'XTickLabel', 2022:2041);
% xlim([1, size(squeeze(incsto), 1)]);
        % Set x-axis ticks and labels for the time range 2022-2035
set(gca, 'XTick', 1:length(years));  % Match the number of years
set(gca, 'XTickLabel', years);  % Update labels to show from 2022 to 2035
xlim([1 length(years)]);  % Explicitly set the x-axis limits to end at 2035

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'FontWeight', 'bold', 'FontSize', 12); 
yline(0.1, '--', 'LineWidth', 1.5, 'Color', 'k'); 
%legend({'Baseline', 'Social protection', 'TPT, recent migrants', 'TPT, pre-entry', ...
%     'TPT, contacts', 'TPT, vulnerables', 'Find and treat, migrants and vulnerables', ...
%     'Find and treat, general population', 'New: TPT+treatment+vaccine'}, 'FontWeight', 'bold', 'FontSize', 12);





incmat   = permute(prctile(incsto,[2.5,50,97.5],2),[2,1,3]);
incmat2  = permute(prctile(incsto2,[2.5,50,97.5],2),[2,1,3]);
incmatRR = permute(prctile(incstoRR,[2.5,50,97.5],2),[2,1,3]);

mrtmat = permute(prctile(mrtsto,[2.5,50,97.5],2),[2,1,3]);


%



ff = figure('Position', [577, 190, 1329, 732]); 
lw = 1.5; 
fs = 14;
% Initialize the plot
hold on;

% Legends for each curve
legendEntries = {'Baseline', 'Social protection', 'TPT, recent migrants', ...
    'TPT, pre-entry', 'TPT, long-term migrants', 'TPT, contacts', ...
    'TPT, vulnerables', 'ACF, migrants and vulnerables', ...
    'ACF, general population', 'New tools: TPT and treatment'};

% Add the zero line, but ensure it's not part of the legend
yline(0.1, '--', 'LineWidth', 1.5, 'Color', 'k', 'HandleVisibility', 'off'); 

% Loop through each scenario to plot the curves one by one
for mi = 1:length(legendEntries)
    plot(squeeze(incsto(:, :, mi)), 'LineWidth', 2); % Plot each scenario
    
    % Update y-axis limits
    yl = ylim; 
    yl(1) = 0; 
    ylim(yl);
    
    % Set x-axis ticks and labels
    set(gca, 'XTick', 1:size(squeeze(incsto), 1));
    set(gca, 'XTickLabel', 2022:2041);
    xlim([1, size(squeeze(incsto), 1)]);

    % Label axes
    xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
    set(gca, 'FontWeight', 'bold', 'FontSize', 12); 
    
    % Update the legend with the currently plotted curves
    legend(legendEntries(1:mi), 'FontWeight', 'bold', 'FontSize', 12);
    
    % Pause to wait for the user to press a key
    pause;
end

hold off;



