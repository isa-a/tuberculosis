clear all; load 15pct;

% Whether to show from 2015, or from 2022
show_extended = 1;

years = 2011:2041;

inc_pct = permute(prctile(incsto, [2.5, 50, 97.5], 2), [2,1,3]);           % Dims: 1.Lo/Md/Hi, 2.Year, 3.Scenario
mrt_pct = permute(prctile(mrtsto, [2.5, 50, 97.5], 2), [2,1,3]);           % Dims: 1.Lo/Md/Hi, 2.Year, 3.Scenario

ff = figure; hold on;
lw = 1.5; tp = 0.1; fs = 14;

cols = linspecer(size(incsto,3));
cols(1, :) = [0, 0, 0.5];
for mi = 1:1
    md   = inc_pct(2,:,mi);
    hilo = inc_pct([1,3],:,mi);

    lg(mi,:) = plot(years, md, 'Color', cols(mi,:), 'LineWidth', lw); hold on;
    jbfill(years, hilo(2,:), hilo(1,:), cols(mi,:), 'None', 1, tp); hold on;

end

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'fontsize', fs);
%ylim([0 10]);
if show_extended
    xlim([years(1) years(end)]);
else
    xlim([2011 2040]);
end
yl = ylim; yl(1) = 0; ylim(yl);

% if show_extended
%     years = [2015:2025];
%     incs1 = [5734,5621,5066,4610,4704,4124,4407,4375,4855,0, 0] * 0.001538;
%     incs1(end-1)   = 9.5;
%     incs1(end)   = 8;
%     plot(years, incs1, '.-', 'MarkerSize',20, 'Color', 0.5*[1 1 1], 'LineWidth',1);
% end

if show_extended
    years = [2011:2025];
    incs1 = [15.59, 15.11, 13.48, 11.91,10.46, 10.17, 9.11, 8.24, 8.36, 7.32, 7.77, 7.65, 8.34, 9.37, 0] / 0.875;
    incs1(end)   = 9.43 / 0.875;
    plot(years, incs1, '.-', 'MarkerSize',20, 'Color', 0.5*[1 1 1], 'LineWidth',1);
end



% if show_extended
%     years = [2015:2025];
%     incs1 = [5734,5621,5066,4610,4704,4124,4407,4375,4855,0, 0] * 0.001538;
%     incs1(end-1)   = 5490*0.001538;
%     incs1(end)   = 5596 * 0.001538;
%     plot(years, incs1, '.-', 'MarkerSize',20, 'Color', 0.5*[1 1 1], 'LineWidth',1);
% end


legend(lg,'Baseline','Improved Tx outcomes', ' + Enhanced TPT, recent migrants','+ Accelerated case-finding (foreign-born)','+ Accelerated case-finding (UK-born)','+ Pre-entry migrant screening, TBI','location','SouthWest');
set(ff,'Position',[344   297   712   563]);


return;

ff = figure; hold on;
lw = 1.5; tp = 0.1; fs = 14;

cols = linspecer(size(incsto,3));
cols(1, :) = [0, 0, 0.5];
for mi = 1:6
    md   = mrt_pct(2,:,mi);
    hilo = mrt_pct([1,3],:,mi);

    lg(mi,:) = plot(years, md, 'Color', cols(mi,:), 'LineWidth', lw); hold on;
    jbfill(years, hilo(2,:), hilo(1,:), cols(mi,:), 'None', 1, tp); hold on;

end

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'fontsize', fs);

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'fontsize', fs);

if show_extended
    xlim([years(1) years(end)]);
else
    xlim([2022 2040]);
end
yl = ylim; yl(1) = 0; ylim(yl);

if show_extended
    years = [2015:2025];
    incs1 = [5734,5621,5066,4610,4704,4124,4407,4375,4855,0, 0] * 0.001538;
    incs1(end-1)   = 5490*0.001538;
    incs1(end)   = 5596 * 0.001538;
    plot(years, incs1, '.-', 'MarkerSize',20, 'Color', 0.5*[1 1 1], 'LineWidth',1);
end

legend(lg,'Baseline','Improved Tx outcomes', ' + Enhanced TPT, recent migrants','+ Accelerated case-finding (foreign-born)','+ Accelerated case-finding (UK-born)','+ Pre-entry migrant screening, TBI','location','SouthWest');
set(ff,'Position',[344   297   712   563]);

return;




% -------------------------------------------------------------------------
% --- Show pie chart of contributions to incidence in final intervention
% --- scenario

tmp1 = mean(incsource,1);
tmp2 = tmp1/sum(tmp1);
% Aggregate all recurrence, still stratified by migrants and domestic
tmp3 = [tmp2(1:4), sum(tmp2([5,7])), sum(tmp2([6,8]))];

figure; pie(fliplr(tmp3),{'','','','','',''});


% -------------------------------------------------------------------------
% --- Plot bar charts showing impact (cases and deaths averted) at each stage
% --- This bit NOT WORKING - need to check

% Find pct cases and deaths averted
allmat  = cat(4, incsto, mrtsto);
callmat = permute(squeeze(sum(allmat(years>=2026,:,:,:),1)),[1,3,2]);      % Dims: 1.Sample 2.Inc/Mor 3.Scenario
pba     = 1 - callmat(:,:,2:end)./callmat(:,:,1);
pba_pct = permute(prctile(pba,[2.5,50,97.5],1),[1,3,2])*100;               % Dims: 1.Lo/Md/Hi 2.Scenario 3.Inc/Mor

lbls = {'Improved Tx outcomes', '+ Enhanced TPT, recent migrants','+ Acc. case-finding \newline (foreign-born)','+ Acc. case-finding \newline (UK-born)','+ Pre-entry migrant screening, TBI'};

tis = {'Incidence','Mortality'};

ff = figure; fs = 14; lw = 1.5;
for ii = 1:2
    subplot(1,2,ii); hold on;
    xx = 1:size(pba_pct,2);
    md = pba_pct(2,:,ii);
    hilo = diff(pba_pct(:,:,ii),[],1);
    
    bar(xx,md);
    errorbar(xx, md, hilo(1,:), hilo(2,:), 'linestyle', 'None', 'linewidth', lw);
    set(gca,'XTick',1:length(xx),'XTickLabel',lbls);
    set(gca,'fontsize',fs);
    title(tis{ii});
end
subplot(1,2,1); ylabel('Percent reduction in cumulative burden, 2026 - 2040');

