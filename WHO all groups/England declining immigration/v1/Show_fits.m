clear all;

% x3=[0.4915,    0.1419 ,   0.2860,    4.3056,    0.9999 ,   0.1770,    6.4304 ,   0.3189  ,  1.0569 ,   1.8553 ,  12.4415 ,   0.9102,    5.5168];


ix0 = size(xsto,1)/2;
nx  = 200;
dx  = round(ix0/nx);
xs  = xsto(ix0:dx:end,:);

% x3=[0.3778,    0.1261  ,  0.0672 ,   4.7561,    1.6836 ,   0.2195,    5.7229 ,   0.3210  ,  1.9110  ,  1.6026 ,  12.8686  ,  0.4666 ,   1.9086];
x3=[0.3778,    0.1261  ,  0.001 ,   1e-4,    1.6836 ,   0.2195,    5.7229 ,   0.3210  ,  1.9110  ,  20 ,  20  ,  0.4666 ,   1.9086];
%xs(ii,:)
for ii = 1:size(x0,1)
    [out, aux] = obj(x0);
    sims(ii,:) = aux.sim;
    inc(:,ii)  = aux.incd;
    pp(ii)     = aux.p_migrect;
end
sim_pct = prctile(sims,[2.5,50,97.5],1);

% Collate data
alldat = [data.incd2010; data.incd2020; data.mort; data.p_migrTB; data.p_migrpopn; data.p_LTBI_inmigr; data.incd_ch2020; data.p_chpopn; data.ch_notifs; data.p_incd_recentinf];
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
set(gca, 'fontsize', 20, 'XTick', 1:size(alldat,1), 'XTickLabel', {'incd2010', 'incd2020', 'mort', 'p_migrTB', 'p_migrpopn', 'p_LTBI', 'incd_ch2020', 'p_chpopn', 'ch_notifs', 'proportion recent'});
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
   histogram(xs(:,ii));
   title(names{ii});
end
