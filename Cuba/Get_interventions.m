clear all; load calibration_res4.mat; load Model_setup.mat;

obj = @(x) get_objective2(x, ref, prm, gps, lhd);

%ix0 = round(size(xsto,1)/2);
%dx  = round(size(xsto,1)/2/150);
%xs  = xsto(ix0:dx:end,:);

[~, maxIndex] = max(outsto);
bestFitParameters = xsto(maxIndex, :);

numSamples = 200;

%  empty  array 2 store dists
distributions = cell(1, length(bestFitParameters));

figure;

% loop thru list of  params
for i = 1:length(bestFitParameters)
 
    centralValue = bestFitParameters(i);
    
    lowerBound = centralValue * 0.5;
    upperBound = centralValue * 1.5;
    
    % generate samples
    samples = lowerBound + (upperBound - lowerBound) * rand(1, numSamples);
    
    distributions{i} = samples;
    
    % plt
    subplot(2, 3, i); 
    histogram(samples, 'Normalization', 'pdf');
    title(sprintf('Parameter %d: Uniform distribution', i));
    xlabel('Value');
    ylabel('Probability Density');
    xlim([lowerBound upperBound]);
end
sgtitle('Uniform Distributions of Best Fit Parameters with ±50% Variation');

% Initialize a matrix to store the sampled values
sampledValues = zeros(numSamples, length(bestFitParameters));

% Take 200 samples in looped order from each parameter's distribution
for j = 1:numSamples
    for i = 1:length(bestFitParameters)
        % Sample one value from each parameter's distribution
        sampledValues(j, i) = distributions{i}(randi(numSamples));
    end
end

opts = odeset('RelTol',1e-9,'AbsTol',1e-9);

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);
    
    init = aux.soln(end,:);

    [p0,r0] = allocate_parameters(xx,p,r,xi);
    M0 = make_model(p0,r0,i,s,gps);

    % ---------------------------------------------------------------------
    % --- Model baseline
    
%     geq = @(t,in) goveqs_basis2(t, in, i, s, M0, agg, sel, r0, p0);
%     [t,soln] = ode15s(geq, [2022:2031], init);
%     sdiff = diff(soln,[],1);
    %incsto(:,ii,1) = sdiff(:,i.aux.inc(1))*1e5;
    
    % ---------------------------------------------------------------------
    % --- Model intervention
    
    p1 = p0; r1 = r0;
    p1.migrTPT = 1;
    M1 = make_model(p1,r1,i,s,gps);
    
    p2 = p0; r2 = r0;
    p2.migrTPT = 1;
    r2.ACF = 0.69*[1 1];
    M2 = make_model(p2,r2,i,s,gps); %acf in every1

    p3 = p0; r3 = r0;
    p3.migrTPT = 1;
    r3.ACF = 0.69*[1 1];
    r3.TPT = 0.69*[0 1];
    M3 = make_model(p3,r3,i,s,gps); % tpt in adults
    
    p4 = p0; r4 = r0;
    p4.migrTPT = 1;
    r4.TPT = 0.69*[1 1];
    r4.ACF = 0.69*[1 1];
    M4 = make_model(p4,r4,i,s,gps); % acf and tpt in every1
    
    models = {M0, M4};
    
    for mi = 1:length(models)
        geq = @(t,in) goveqs_scaleup(t, in, i, M0, models{mi}, [2024 2029], agg, sel, r, p0);
        [t,soln] = ode15s(geq, [2022:2041], init, opts);
        
        sdiff = diff(soln,[],1);
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1))*1e5;
        
        % Get proportions from different sources
        vec = sdiff(end,i.aux.incsources)*1e5;  
        props(ii,:,mi) = vec/sum(vec);
    end
end
fprintf('\n');













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





opts = odeset('RelTol',1e-9,'AbsTol',1e-9);

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    
    if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
    
    xx = xs(ii,:);
    [out,aux] = obj(xx);
    
    init = aux.soln(end,:);

    [p0,r0] = allocate_parameters(xx,p,r,xi);
    M0 = make_model(p0,r0,i,s,gps);

    % ---------------------------------------------------------------------
    % --- Model baseline
    
%     geq = @(t,in) goveqs_basis2(t, in, i, s, M0, agg, sel, r0, p0);
%     [t,soln] = ode15s(geq, [2022:2031], init);
%     sdiff = diff(soln,[],1);
    %incsto(:,ii,1) = sdiff(:,i.aux.inc(1))*1e5;
    
    % ---------------------------------------------------------------------
    % --- Model intervention
    
    p1 = p0; r1 = r0;
    p1.migrTPT = 1;
    M1 = make_model(p1,r1,i,s,gps);
    
    p2 = p0; r2 = r0;
    p2.migrTPT = 1;
    r2.ACF = 0.69*[1 1];
    M2 = make_model(p2,r2,i,s,gps); %acf in every1

    p3 = p0; r3 = r0;
    p3.migrTPT = 1;
    r3.ACF = 0.69*[1 1];
    r3.TPT = 0.69*[0 1];
    M3 = make_model(p3,r3,i,s,gps); % tpt in adults
    
    p4 = p0; r4 = r0;
    p4.migrTPT = 1;
    r4.TPT = 0.69*[1 1];
    r4.ACF = 0.69*[1 1];
    M4 = make_model(p4,r4,i,s,gps); % acf and tpt in every1
    
    models = {M0, M4};
    
    for mi = 1:length(models)
        geq = @(t,in) goveqs_scaleup(t, in, i, M0, models{mi}, [2024 2029], agg, sel, r, p0);
        [t,soln] = ode15s(geq, [2022:2041], init, opts);
        
        sdiff = diff(soln,[],1);
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1))*1e5;
        
        % Get proportions from different sources
        vec = sdiff(end,i.aux.incsources)*1e5;  
        props(ii,:,mi) = vec/sum(vec);
    end
end
fprintf('\n');


mat = permute(prctile(incsto,[2.5,50,97.5],2),[2,1,3]);

cols = linspecer(size(mat,3));
figure; lw = 3; fs = 14;

xx = [2022:2040];
for ii = 1:size(mat,3)
   plt = mat(:,:,ii);
   lg(ii,:) = plot(xx, plt(2,:), 'Color', cols(ii,:), 'linewidth', lw); hold on;
   jbfill(xx, plt(3,:), plt(1,:), cols(ii,:), 'None', 1, 0.1); hold on;
end
yl = ylim; yl(1) = 0; ylim(yl);
set(gca,'fontsize',fs);
ylabel('Incidence per 100,000 population');
xlabel('Year');
yline(0.1,'k--','LineWidth', 2);
yline(1,'k--','LineWidth', 2);
title('Sampled params');

legend(lg, 'Baseline','ACF in everyone','ACF + adults TPT','ACF + TPT in everyone','location','SouthWest');


return;
%




paramset = xsto(6209, :);



[out, aux] = obj(paramset);
init = aux.soln(end, :);

[p0, r0] = allocate_parameters(paramset, p, r, xi);
M0 = make_model(p0, r0, i, s, gps);

% ---------------------------------------------------------------------
% --- Model baseline
% geq = @(t, in) goveqs_basis2(t, in, i, s, M0, agg, sel, r0, p0);
% [t, soln] = ode15s(geq, [2022:2031], init);
% sdiff = diff(soln, [], 1);
% incsto(:, 1, 1) = sdiff(:, i.aux.inc(1)) * 1e5;

% ---------------------------------------------------------------------
% --- Model intervention

p1 = p0; r1 = r0;
p1.migrTPT = 1;
M1 = make_model(p1, r1, i, s, gps);

p2 = p0; r2 = r0;
p2.migrTPT = 1;
r2.ACF = 0.69 * [1 1];
M2 = make_model(p2, r2, i, s, gps); % ACF in everyone

p3 = p0; r3 = r0;
p3.migrTPT = 1;
r3.ACF = 0.69 * [1 1];
r3.TPT = 0.69 * [0 1];
M3 = make_model(p3, r3, i, s, gps); % TPT in adults

p4 = p0; r4 = r0;
p4.migrTPT = 1;
r4.TPT = 0.69 * [1 1];
r4.ACF = 0.69 * [1 1];
M4 = make_model(p4, r4, i, s, gps); % ACF and TPT in everyone

models = {M0, M4};

opts = odeset('RelTol',1e-14,'AbsTol',1e-14);

incsto = zeros(2041 - 2022, 1, length(models)); 
props = zeros(1, length(i.aux.incsources), length(models));

for mi = 1:length(models)
    geq = @(t, in) goveqs_scaleup(t, in, i, M0, models{mi}, [2024 2029], agg, sel, r, p0);
    [t, soln] = ode15s(geq, [2022:2041], init, opts);

    sdiff = diff(soln, [], 1);
    incsto(:, 1, mi) = sdiff(:, i.aux.inc(1)) * 1e5;


    vec = sdiff(end, i.aux.incsources) * 1e5;
    props(1, :, mi) = vec / sum(vec);
end

cols = linspecer(length(models));
figure; lw = 3; fs = 14;

xx = 2022:2040;
for ii = 1:length(models)
    plt = incsto(:, 1, ii);
    lg(ii, :) = plot(xx, plt, 'Color', cols(ii, :), 'linewidth', lw); hold on;
end

yl = ylim; yl(1) = 0; ylim(yl);
set(gca, 'fontsize', fs);
ylabel('Incidence per 100,000 population');
xlabel('Year');
yline(0.1, 'k--', 'LineWidth', 2);
yline(1, 'k--', 'LineWidth', 2);

legend(lg, 'Baseline', 'ACF in everyone', 'ACF + adults TPT', 'ACF + TPT in everyone', 'location', 'SouthWest');
title('Single parameter set');




%




