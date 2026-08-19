clear all; load intvn_res2026_04;
years = 2015:2041;
% --- Cumulative cases/deaths averted vs baseline, from 2026 onwards ---
allmat  = cat(4, incsto, mrtsto);
% Sum cumulative burden from 2026, dims: 1.Sample 2.Inc/Mor 3.Scenario
callmat = permute(squeeze(sum(allmat(years >= 2026, :, :, :), 1)), [1, 3, 2]);
% Percent reduction vs baseline (scenario 1)
pba     = 1 - callmat(:, :, 2:end) ./ callmat(:, :, 1);
pba_pct = permute(prctile(pba, [2.5, 50, 97.5], 1), [1, 3, 2]) * 100;
% Dims: 1.Lo/Md/Hi  2.Scenario(2-6 vs baseline)  3.Inc/Mor
lbls = {
    'Improved Tx outcomes', ...
    '+ Enhanced TPT (migrants)', ...
    '+ Case-finding (foreign-born)', ...
    '+ Case-finding (UK-born)', ...
    '+ Pre-entry screening (TBI)'
};
tis  = {'Incidence', 'Mortality'};
cols = linspecer(5);
ff = figure; fs = 13; lw = 1.5;
set(ff, 'Position', [200 200 1000 480]);
for ii = 1:2
    subplot(1, 2, ii); hold on;
    xx   = 1:size(pba_pct, 2);
    md   = pba_pct(2, :, ii);
    lo   = md - pba_pct(1, :, ii);
    hi   = pba_pct(3, :, ii) - md;
    % Coloured bars
    for bi = 1:length(xx)
        bar(xx(bi), md(bi), 0.6, 'FaceColor', cols(bi,:), 'EdgeColor', 'k');
    end
    % 95% CI whiskers
    errorbar(xx, md, lo, hi, 'linestyle', 'None', 'linewidth', lw, 'Color', 'k', 'CapSize', 8);
    set(gca, 'XTick', xx, 'XTickLabel', lbls, ...
        'FontSize', fs, 'XTickLabelRotation', 25, 'FontWeight', 'bold');
    title(tis{ii}, 'FontSize', fs+1, 'FontWeight', 'bold');
    grid on; box on;
    yl = ylim; yl(1) = 0; ylim(yl);
end
subplot(1, 2, 1);
ylabel('% reduction in cumulative burden (2026–2041)', 'FontSize', fs, 'FontWeight', 'bold');
subplot(1, 2, 2);
ylabel('% reduction in cumulative burden (2026–2041)', 'FontSize', fs, 'FontWeight', 'bold');