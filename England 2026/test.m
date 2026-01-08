clear all;

load 9pct
incsto_A = incsto;
mrtsto_A = mrtsto;

load 15pct
incsto_B = incsto;
mrtsto_B = mrtsto;


show_extended = 1;

years = 2011:2041;

inc_pct_A = permute(prctile(incsto_A, [2.5, 50, 97.5], 2), [2,1,3]);
mrt_pct_A = permute(prctile(mrtsto_A, [2.5, 50, 97.5], 2), [2,1,3]);

inc_pct_B = permute(prctile(incsto_B, [2.5, 50, 97.5], 2), [2,1,3]);
mrt_pct_B = permute(prctile(mrtsto_B, [2.5, 50, 97.5], 2), [2,1,3]);


nyA = size(inc_pct_A, 2);
nyB = size(inc_pct_B, 2);
ny  = min([nyA, nyB, length(years)-1]);   
years_inc = years(1:ny);

ff = figure; hold on;
lw = 1.5; tp = 0.1; fs = 14;


cols = [0 0 0.5];
for mi = 1:1
    md   = inc_pct_A(2,1:ny,mi);
    hilo = inc_pct_A([1,3],1:ny,mi);

    lg(1,:) = plot(years_inc, md, 'Color', cols, 'LineWidth', lw); hold on;
    jbfill(years_inc, hilo(2,:), hilo(1,:), cols, 'None', 1, tp); hold on;
end


cols = [0.8 0 0];
for mi = 1:1
    md   = inc_pct_B(2,1:ny,mi);
    hilo = inc_pct_B([1,3],1:ny,mi);

    lg(2,:) = plot(years_inc, md, 'Color', cols, 'LineWidth', lw); hold on;
    jbfill(years_inc, hilo(2,:), hilo(1,:), cols, 'None', 1, tp); hold on;
end

xlabel('Year', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Rate per 100,000 population', 'FontWeight', 'bold', 'FontSize', 12);
set(gca, 'fontsize', fs);

%ylim([0 10]);
if show_extended
    xlim([2011 years(end)]);
else
    xlim([2022 2040]);
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
    incs1(end) = 9.43 / 0.875;
    plot(years, incs1, '.-', 'MarkerSize',20, ...
        'Color', 0.5*[1 1 1], 'LineWidth',1);

    % --- FIX 2: restore years (no rename, just put it back) ---
    years = 2015:2041;
end

legend(lg, 'Constant LTBI positivity throughout', 'Increased LTBI post-2020', ...
       'location', 'SouthWest');

set(ff,'Position',[344 297 712 563]);

return;
