clear all; load model_fits2; load Model_setup.mat;
obj = @(x) get_objective(x, prm, ref, sel, agg, gps, lhd);

ix0 = 2e4; nx = 250; dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);

% inds = find(outsto==max(outsto));
% xs = xsto(inds(1),:);

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 
    [out, aux] = obj(xs(ii,:));
    sims(ii,:) = aux.sim;
    inct(:,:,ii) = diff(aux.soln(:,i.aux.inc),1)*1e5;
end
fprintf('\n');
sim_pct = prctile(sims,[2.5,50,97.5],1);
inc_pct = permute(prctile(inct,[2.5,50,97.5],3),[3,1,2]);                  % 1. Lo/Md/Hi, 2.Time, 3.All/HIV only

alldat  = [data.inc2019_all; data.inc2019_H1; data.notifs(1)*[0.85 1 1.15]; data.ART_covg; data.HIV_prev; data.mort2019_H0; data.mort2019_H1; data.psym]';
den = alldat(2,:);

%inc_all, inc_h1, noti, ART_covg, HIV_prev, mort_H0, mort_H1, psym

% Normalise all with respect to central data estimates
plt_sim = sim_pct./den;
plt_dat = alldat./den;


% --- Show all on a plot --------------------------------------------------

figure; fs = 14;

% --- Incidence plot
% subplot(1,2,1);

plt = sum(inc_pct,3);
plot(plt(2,:)); hold on;
jbfill(1:size(plt,2),plt(3,:),plt(1,:),'b','None',1,0.3); hold on;
yl = ylim; yl(1) = 0; ylim(yl);

plt = sum(inc_pct(:,:,[2,3]),3);
plot(plt(2,:)); hold on;
jbfill(1:size(plt,2),plt(3,:),plt(1,:),'r','None',1,0.3); hold on;
yl = ylim; yl(1) = 0; ylim(yl);

ylabel('Incidence');
set(gca,'fontsize',fs);

% --- Comparing model outputs with simulations
figure; hold on;

allplt = cat(3,plt_dat,plt_sim);
cols = {'r','b'};

for ii = 1:2
    plt = allplt(:,:,ii);
    hilo = diff(plt,1); md = plt(2,:);
    xpts = [1:length(plt)] + (-1)^ii*0.1;
    plot(xpts, md, '.', 'Color', cols{ii}, 'markersize', 24);
    errorbar(xpts, md, hilo(1,:), hilo(2,:), 'LineStyle', 'None', 'Color', cols{ii});
end
set(gca,'fontsize',fs,'XTick',1:size(sims,2),'XTickLabel',{'Incd all', 'Incd H1', 'Noti', 'ART covg', 'HIV prev', 'Mort H0', 'Mort H1', 'Symp'});
xtickangle(45);

yl = ylim; yl(1) = 0; ylim(yl);

% --- Show posterior densities
figure;
names  = {}; fnames = fieldnames(xi);
for ii = 1:length(fnames)
    fname = fnames{ii};
    names = [names, repmat({fname},1,length(xi.(fname)))];
end
for ii = 1:12
   subplot(4,3,ii); 
   histogram(xsto(:,ii));
   title(names{ii});
end

% --- Show cross-correlations
figure; lim = 1;
for ir = 1:12
    for ic = 1:12
        subplot(12,12,lim);
        plot(xs(:,ir), xs(:,ic), '.');
        if ic==1
            ylabel(names{ir});
        end
        if ir==12
            xlabel(names{ic});
        end
        lim=lim+1;
    end
end

