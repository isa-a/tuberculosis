clear all; load intvn_res;

years = 2022:2041;

inc_pct = permute(prctile(incsto, [2.5, 50, 97.5], 2), [2,1,3]);           % Dims: 1.Lo/Md/Hi, 2.Year, 3.Scenario
inc_pct = inc_pct(:,:,[1,3:6,2]);

mrt_pct = permute(prctile(mrtsto, [2.5, 50, 97.5], 2), [2,1,3]);           % Dims: 1.Lo/Md/Hi, 2.Year, 3.Scenario
mrt_pct = mrt_pct(:,:,[1,3:6,2]);

ff = figure; hold on;
lw = 1.5; tp = 0.1; fs = 14;

cols = linspecer(size(incsto,3));
cols(1, :) = [0, 0, 0.5];
for mi = 1:6
    md   = inc_pct(2,:,mi);
    hilo = inc_pct([1,3],:,mi);

    lg(mi,:) = plot(years, md, 'Color', cols(mi,:), 'LineWidth', lw); hold on;
    jbfill(years, hilo(2,:), hilo(1,:), cols(mi,:), 'None', 1, tp); hold on;

end

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
% ylim([0, 11]);
set(gca, 'fontsize', fs);
xlim([years(1) years(end)]);
ylim([0 12]);
%yl = ylim; yl(1) = 0; ylim(yl);


legend(lg,'Baseline','Enhanced TPT, recent migrants','+ Improved Tx outcomes','+ Accelerated case-finding (foreign-born)','+ Accelerated case-finding (UK-born)','+ Pre-entry migrant screening, TBI');

set(ff,'Position',[577, 190, 1029, 732]);

return;
%

gg = figure; hold off;
lw = 1.5; tp = 0.1; fs = 14;

cols2 = linspecer(size(mrtsto,3));


for mi = 1:6
    md   = mrt_pct(2,:,mi);
    hilo = mrt_pct([1,3],:,mi);

    lg(mi,:) = plot(years, md, 'Color', cols2(mi,:), 'LineWidth', lw); hold on;
    jbfill(years, hilo(2,:), hilo(1,:), cols2(mi,:), 'None', 1, tp); hold on;

end

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
% ylim([0, 11]);
set(gca, 'fontsize', fs);
xlim([years(1) years(end)]);

yl = ylim; yl(1) = 0; ylim(yl);


legend(lg,'Baseline','Enhanced TPT, recent migrants','+ Improved Tx outcomes','+ Accelerated case-finding (foreign-born)','+ Accelerated case-finding (UK-born)','+ Pre-entry migrant screening, TBI');

set(gg,'Position',[577, 190, 1029, 732]);