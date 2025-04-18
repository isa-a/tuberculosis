clear all;

% x3=[0.4915,    0.1419 ,   0.2860,    4.3056,    0.9999 ,   0.1770,    6.4304 ,   0.3189  ,  1.0569 ,   1.8553 ,  12.4415 ,   0.9102,    5.5168];
load optim_res3;

x3(xi.betadec)  = 0;
% x3(xi.gamma(2)) = 2;

ix0 = size(xsto,1)/2;
nx  = 200;
dx  = round(ix0/nx);
xs  = xsto(ix0:dx:end,:);

%xs(ii,:)
for ii = 1:size(x3,1)
    [out, aux, msg] = obj(x3);
    sims(ii,:) = aux.sim;
    inc(:,ii)  = aux.incd;
    pp(ii)     = aux.p_migrect;
end
sim_pct = prctile(sims,[2.5,50,97.5],1);

% aux.p_incd_recentinf

% Collate data
alldat = [data.incd2010; data.incd2020; data.mort; data.p_migrTB; data.p_migrpopn; data.p_LTBI_inmigr; data.p_incd_recentinf];
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
set(gca, 'fontsize', 20, 'XTick', 1:size(alldat,1), 'XTickLabel', {'incd2010', 'incd2020', 'mort', 'p_migrTB', 'p_migrpopn', 'p_LTBI', 'proportion recent'});
yl = ylim; yl(1) = 0; ylim(yl);


% Plot incidence timeseries
figure; hold on;
plot(2010:2029, inc,'Color',0.5*[1 1 1]);
yl = ylim; yl(1) = 0; ylim(yl);
xlim([2010, 2024])

return;

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
   histogram(xs(:,ii));
   title(names{ii});
end
