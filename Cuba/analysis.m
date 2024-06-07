clear all; load calibration_res5.mat; load Model_setup.mat;

obj = @(x) get_objective2(x, ref, prm, gps, lhd);
ix0 = round(size(xsto,1)/2);
dx  = round(size(xsto,1)/2/200);
xs  = xsto(ix0:dx:end,:);

opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

timereached = NaN(size(xs, 1), 1); % empty for time reached
ch_inc = NaN(size(xs, 1), 1); % 


mk = round(size(xs, 1) / 25);
for ii = 1:size(xs, 1)
    
    if mod(ii, mk) == 0
        fprintf('%0.5g ', ii / mk); 
    end
    
    xx = xs(ii, :);
    [out, aux] = obj(xx);
    
    init = aux.soln(end, :);

    [p0, r0] = allocate_parameters(xx, p, r, xi);
    M0 = make_model(p0, r0, i, s, gps);

    p4 = p0; r4 = r0;
    p4.migrTPT = 1;
    r4.TPT = 3.91202 * [1 1];
    r4.ACF = 3.91202 * [1 1];
    M4 = make_model(p4, r4, i, s, gps); % ACF and TPT in everyone
    
    models = {M0, M4};
    
    for mi = 1:length(models)
        geq = @(t, in) goveqs_scaleup(t, in, i, M0, models{mi}, [2024 2029], agg, sel, r, p0);
        [t, soln] = ode15s(geq, [2022:2100], init, opts);
        
        sdiff = diff(soln, [], 1);
        
        incsto(:,ii,mi) = sdiff(:,i.aux.inc(1));
        incsto2(:,ii,mi) = sdiff(:,i.aux.inc(2));
        incsto3(:,ii,mi) = sdiff(:,i.aux.inc(3));

        % find  year when inc reaches 1 per mill
        idx = find(sdiff(:, i.aux.inc(1)) * 1e5 <= 0.1, 1); % 1 per million = 0.1 per 100,000
        if ~isempty(idx)
            timereached(ii, mi) = t(idx);
            if mi == 2 %  model 4 is second model in the list
                ch_inc(ii, 1) = soln(idx, i.aux.inc(2)); % *1e5 for per 100,000 children
            end

        else
            timereached(ii, mi) = NaN; % Set to NaN if threshold is never reached
        end
    end
end

fprintf('\n');

samplescell = mat2cell(xs, ones(size(xs, 1), 1), size(xs, 2));
resultstable = table(samplescell, timereached, ch_inc, 'VariableNames', {'Parameter Set', 'Year reached', 'Child incidence at this point in time'});
disp(resultstable);

incsto = incsto*1e5;
incsto2 = incsto2*1e5;
incsto3 = incsto3*1e5;

%parameter set that gives latest elimination
figure;
subplot(1,3,1);
plot(2022:2099, squeeze(incsto2(:,88,:)))
yline(0.1, 'k--', 'LineWidth', 2); %children

subplot(1,3,2);
plot(2022:2099, squeeze(incsto3(:,88,:)))
yline(0.1, 'k--', 'LineWidth', 2); %adults

subplot(1,3,3);
plot(2022:2099, squeeze(incsto(:,88,:)))
yline(0.1, 'k--', 'LineWidth', 2); % all incidence



%  1st and last year of elimination
[earliest_yr, earliest_idx] = min(timereached(:, 2));
[latest_yr, latest_idx] = max(timereached(:, 2));
% params
earliest_params = xs(earliest_idx, :);
latest_params = xs(latest_idx, :);

disp('Earliest Year Parameter Set:');
disp(earliest_params);
disp(['Earliest Year: ', num2str(earliest_yr)]);

disp('Latest Year Parameter Set:');
disp(latest_params);
disp(['Latest Year: ', num2str(latest_yr)]);

% compare
comp_tab = table(earliest_params', latest_params', 'VariableNames', {'Earliest Year Params', 'Latest Year Params'});
disp('Comparison of Parameter Sets:');
disp(comp_tab);


% Create the matrix with parameter sets and child incidence
parameterMatrix = [xs, ch_inc];

% Perform partial correlation analysis
% Removing rows with NaN values
parameterMatrix = parameterMatrix(~any(isnan(parameterMatrix), 2), :);
[partialCorrMatrix, pValues] = partialcorr(parameterMatrix);

fprintf('\n');

% Display results
disp('Partial Correlation Matrix:');
disp(partialCorrMatrix);


% Optionally, visualize the partial correlation matrix
figure;
imagesc(partialCorrMatrix);
colorbar;
title('Partial Correlation Matrix');
xlabel('Parameters and Child Incidence');
ylabel('Parameters and Child Incidence');