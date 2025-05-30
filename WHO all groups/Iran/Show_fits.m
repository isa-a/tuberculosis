clear all; %load optim_resUK3.mat;
% load optim_res_noVULN5_v2.mat;
load optim_res_MAIN_OMAN;



ix0 = size(xsto,1)/2;
nx  = 200;
dx  = round(ix0/nx);
xs  = xsto(ix0:dx:end,:);


for ii = 1:size(xs,1)
    [out, aux] = obj(xs(ii,:));
    sims(ii,:) = aux.sim;
    inc(:,ii)  = aux.incd;
    pp(ii)     = aux.p_migrect;
end
sim_pct = prctile(sims,[2.5,50,97.5],1);

% Collate data
%alldat = [data.incd2010; data.incd2020; data.mort; data.p_migrTB; data.p_migrpopn; data.p_LTBI_inmigr; data.propincd_ch; data.p_chpopn; data.p_adpopn; data.ch_notifs; data.vuln_prev; data.vuln_relrisk;];
alldat = [data.incd2010; data.incd2020; data.mort; data.p_migrTB; data.p_migrpopn; data.p_LTBI_inmigr; data.p_vulnpopn; data.p_vulnTB; data.propincd_ch; data.p_chpopn; data.p_adpopn; data.ch_notifs];
den = alldat(:,2)';

sim_plt = sim_pct./den;
dat_plt = alldat'./den;

figure; ms = 24; hold on;

md = dat_plt(2,:); hilo = diff(dat_plt,[],1);
xx = (1:length(md)) - 0.1;
plot(xx, md, 'r.', 'markersize',ms);
errorbar(xx, md, hilo(1,:), hilo(2,:), 'Color', 'r', 'linestyle', 'None');

md = sim_plt(2,:); hilo = diff(sim_plt,[],1);
xx = (1:length(md)) + 0.1;
plot(xx, md, 'b.', 'markersize',ms);
errorbar(xx, md, hilo(1,:), hilo(2,:), 'Color', 'b', 'linestyle', 'None');
%set(gca,'XTick',1:size(alldat,1),'XTickLabel',fieldnames(data));
set(gca, 'fontsize', 20, 'XTick', 1:size(alldat,1), 'XTickLabel', {'incd2010', 'incd2020', 'mort', 'p_migrTB', 'p_migrpopn', 'p_LTBI',  'p_vulnpopn', 'p_vulnTB', 'incd_ch2020', 'p_chpopn', 'p_adpopn', 'ch_notifs'});
yl = ylim; yl(1) = 0; ylim(yl);

return

% Plot incidence timeseries
figure; hold on;
plot(inc,'Color',0.5*[1 1 1]);



% Plot data for comparison
vecs = [data.incd2010; data.incd2020]';
hilo = diff(vecs,[],1);
plot(1, vecs(2,1),'r.','markersize',24);
plot(10,vecs(2,2),'r.','markersize',24);
errorbar([1, 10], vecs(2,:), hilo(1,:), hilo(2,:),'linestyle','None');
yl = ylim; yl(1) = 0; ylim(yl);

% Plot posterior densities
figure;
names  = {};
fnames = fieldnames(xi);
for ii = 1:length(fnames)
    fname = fnames{ii};
    names = [names, repmat({fname},1,length(xi.(fname)))];
end
for ii = 1:16
   subplot(4,4,ii); 
   histogram(x0(:,ii));
   title(names{ii});
end


for ii = 1:size(x0sto,1)
    [out, aux] = obj(x0sto(ii,:));
    sims(ii,:) = aux.sim;
    inc(:,ii)  = aux.incd;
    pp(ii)     = aux.p_migrect;
end

% Collate data
alldat = [data.incd2010; data.incd2020; data.mort; data.p_migrTB; data.p_migrpopn; data.p_LTBI_inmigr; data.p_vulnpopn; data.p_vulnTB; data.propincd_ch; data.p_chpopn; data.p_adpopn; data.ch_notifs];
den = alldat(:,2)';

% Compute the simulation results for plotting
sim_plt = sims ./ den;
dat_plt = alldat' ./ den;

% Plotting
figure; ms = 24; hold on;

% Plot the real data with red error bars
md = dat_plt(2,:);
hilo = diff(dat_plt, [], 1);
xx = (1:length(md)) - 0.1;
plot(xx, md, 'r.', 'markersize', ms);
errorbar(xx, md, hilo(1,:), hilo(2,:), 'Color', 'r', 'linestyle', 'None');

% Define colors for each parameter set
colors = {'b', 'k', 'g', 'y', 'm'}; % Blue, black, green

% Plot simulation results for each parameter set
for ii = 1:size(sim_plt, 1)
    md = sim_plt(ii, :);
    xx = (1:length(md)) + 0.1 + 0.1 * (ii - 1); % Adjust x positions slightly to avoid overlap
    plot(xx, md, '.', 'Color', colors{ii}, 'markersize', ms);
end

% Set X-axis labels
set(gca, 'fontsize', 20, 'XTick', 1:size(alldat,1), 'XTickLabel', {'incd2010', 'incd2020', 'mort', 'p_migrTB', 'p_migrpopn', 'p_LTBI', 'incd_ch2020', 'p_chpopn', 'p_adpopn', 'ch_notifs',  'vuln_prev', 'vuln_relrisk'});
yl = ylim; yl(1) = 0; ylim(yl);

% Create legend
h1 = plot(nan, nan, 'r.', 'markersize', ms);  % Placeholder for real data
h2 = plot(nan, nan, 'b.', 'markersize', ms);  % Placeholder for Set 1
h3 = plot(nan, nan, 'k.', 'markersize', ms);  % Placeholder for Set 2
h4 = plot(nan, nan, 'g.', 'markersize', ms);  % Placeholder for Set 3
h5 = plot(nan, nan, 'y.', 'markersize', ms);  % Placeholder for Set 2
h6 = plot(nan, nan, 'm.', 'markersize', ms);  % Placeholder for Set 3
legend([h1, h2, h3, h4], {'Real Data', 'Set 1', 'Set 2', 'Set 3'}, 'Location', 'Best');




for ii = 1:size(x0sto,1)
    [out, aux] = obj(x0sto(ii,:));
    sims(ii,:) = aux.sim;
    inc(:,ii)  = aux.incd;
    pp(ii)     = aux.p_migrect;
end

% Collate data
alldat = [data.incd2010; data.incd2020; data.mort; data.p_migrTB; data.p_migrpopn; data.p_LTBI_inmigr; data.p_vulnpopn; data.p_vulnTB; data.propincd_ch; data.p_chpopn; data.p_adpopn; data.ch_notifs];
den = alldat(:,2)';

% Compute the simulation results for plotting
sim_plt = sims ./ den;
dat_plt = alldat' ./ den;

% Plotting
figure; ms = 24; hold on;

% Plot the real data with red error bars
md = dat_plt(2,:);
hilo = diff(dat_plt, [], 1);
xx = (1:length(md)) - 0.1;
plot(xx, md, 'r.', 'markersize', ms);
errorbar(xx, md, hilo(1,:), hilo(2,:), 'Color', 'r', 'linestyle', 'None');

% Define colors for each parameter set (using a larger set of colors for 10 sets)
colors = lines(10); % Automatically selects distinct colors for up to 10 plots

% Plot simulation results for each parameter set
for ii = 1:size(sim_plt, 1)
    md = sim_plt(ii, :);
    xx = (1:length(md)) + 0.1 + 0.1 * (ii - 1); % Adjust x positions slightly to avoid overlap
    plot(xx, md, '.', 'Color', colors(ii,:), 'markersize', ms);
end

% Set X-axis labels
set(gca, 'fontsize', 20, 'XTick', 1:size(alldat,1), 'XTickLabel', {'incd2010', 'incd2020', 'mort', 'p_migrTB', 'p_migrpopn', 'p_LTBI',  'p_vulnpopn', 'p_vulnTB', 'incd_ch2020', 'p_chpopn', 'p_adpopn', 'ch_notifs'});
yl = ylim; yl(1) = 0; ylim(yl);

% Create legend (automatically adapts for 10 sets)
h1 = plot(nan, nan, 'r.', 'markersize', ms);  % Placeholder for real data
legend_labels = {'Real Data'};
for ii = 1:10
    h(ii) = plot(nan, nan, '.', 'Color', colors(ii,:), 'markersize', ms); % Placeholder for each simulation
    legend_labels{ii+1} = ['Set ' num2str(ii)];
end
legend([h1, h], legend_labels, 'Location', 'Best');

