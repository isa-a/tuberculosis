clear all; load optim_res_MAIN.mat; load Model_setup;

obj = @(x) get_objective2(x, ref, prm, gps, prm.contmat, lhd);

set1 = {'ds'};
set2 = {'dom','migr','vuln'};
set3 = {'L','P','R','T'};
[inci, incs, incd, lim] = get_addresses({set3, set2, set1}, [], [], [], 0);
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

midpt = true; 
if midpt
    xs = x0sto(5,:);
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

    % ---------------------------------------------------------------------
    % --- Combined Model Intervention
    
    TPTcov = -log(1-0.99); 
    ACFcov = -log(1-0.99); 

    % Apply all interventions directly to r0 and p0 for a combined effect
    r0.TPT = TPTcov * [1 1 1 1];       % Apply TPT to all target groups
    r0.progression = r0.progression * 0.9; % For TPT in contacts
    r0.reactivation = r0.reactivation * 0.9;
    r0.ACF = ACFcov * [1 1 1 1];       % Apply ACF to all target groups
    r0.ACF2 = r0.ACF;
    p0.TPTeff = 0.8;                   % Set TPT effectiveness
    p0.migrTPT = 1;                    % TPT for migrants at point of entry

    M_combined = make_model(p0, r0, i, s, gps, prm.contmat);

    models = {M0, M_combined};  % Baseline and combined intervention model
    
    for mi = 1:length(models)
        
        geq = @(t,in) goveqs_scaleup(t, in, i, s, M0, models{mi}, p0, p0, [2024 2029], agg, sel, r0);
        [t,soln] = ode15s(geq, [2022:2041], init, opts);
        
        endsolsto(mi,:) = soln(end,:);
        
        sdiff = diff(soln,[],1);
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1))*1e5;
        mrtsto(:,ii,mi) = sdiff(:,i.aux.mort)*1e5;
        
        vec = sdiff(end,i.aux.incsources)*1e5;
        props(ii,:,mi) = vec/sum(vec);
        vec1 = vec;
    end
end
fprintf('\n');

% Plotting baseline and combined intervention incidence rates

ff = figure('Position', [577, 190, 1029, 732]); 
lw = 1.5; 
fs = 14;
hold on;

% Define legend entries
legendEntries = {'Baseline', 'Combined Interventions'};

% Define the years for the x-axis, from 2022 to 2035
years = 2022:2035;

% Plot the baseline incidence (M0)
plot(years, mean(incsto(1:length(years), :, 1), 2), 'LineWidth', lw, 'DisplayName', legendEntries{1});

% Plot the combined intervention incidence
plot(years, mean(incsto(1:length(years), :, 2), 2), 'r-', 'LineWidth', lw, 'DisplayName', legendEntries{2});

% Formatting
xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'FontWeight', 'bold', 'FontSize', 12); 
xlim([years(1) years(end)]);
ylim([0 max(max(mean(incsto(1:length(years), :, 1), 2)), max(mean(incsto(1:length(years), :, 2), 2))) * 1.1]);

% Legend
legend('show', 'Location', 'best', 'FontWeight', 'bold', 'FontSize', 12);

hold off;
