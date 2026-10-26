clear all; load intvn_res2026_04;
years = 2015:2041;

% Find pct cases and deaths averted
allmat  = cat(4, incsto, mrtsto);
callmat = permute(squeeze(sum(allmat(years >= 2026, :, :, :), 1)), [1, 3, 2]);
pba     = 1 - callmat(:, :, 2:end) ./ callmat(:, :, 1);
pba_pct = permute(prctile(pba, [2.5, 50, 97.5], 1), [1, 3, 2]) * 100;


lbls = {
'Improved treatment outcomes', ...
'+ Enhanced TPT (migrants)', ...
'+ Case-finding (migrants)', ...
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
    hilo = diff(pba_pct(:, :, ii), [], 1);


    for bi = 1:length(xx)
        bar(xx(bi), md(bi), 0.6, 'FaceColor', cols(bi,:), 'EdgeColor', 'k');
    end

    errorbar(xx, md, hilo(1,:), hilo(2,:), 'linestyle', 'None', 'linewidth', lw, 'Color', 'k', 'CapSize', 8);
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


disp(array2table(squeeze(pba_pct(:,:,1))', ...
    'VariableNames', {'Lo','Mid','Hi'}, 'RowNames', lbls))