% Load data
data = readtable('observed_ratios_with_incidence.xlsx');

% Get unique countries (sorted alphabetically to match previous plot)
countries = unique(data.Country);
n_countries = length(countries);

% Create figure
figure('Position', [100, 100, 1200, 600]);
hold on;

% Generate same colors as previous plot (must use same logic)
if n_countries <= 12
    colors = lines(n_countries);
else
    % Combine multiple colormaps for better distinction
    colors1 = lines(12);           % First 12 colors
    colors2 = hsv(12);             % Next 12 colors
    colors3 = jet(n_countries-24); % Remaining colors
    colors = [colors1; colors2; colors3];
end

% Plot each country
for i = 1:n_countries
    % Filter data for this country
    idx = strcmp(data.Country, countries{i});
    country_data = data(idx, :);
    
    % Plot with specific color (lines connecting points from same country)
    plot(country_data.Incidence_per_100k, country_data.Ratio, '-o', ...
         'LineWidth', 2, ...
         'MarkerSize', 6, ...
         'Color', colors(i,:), ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:), ...
         'DisplayName', countries{i});
end

% Add horizontal line at 1.37 (median)
yline(1.37, '--r', 'LineWidth', 2, 'DisplayName', 'Median Ratio = 1.37');

% Labels and formatting
xlabel('Incidence per 100k', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Scaling factor', 'FontSize', 12, 'FontWeight', 'bold');
% title('LTBI Scaling Ratios vs TB Incidence by Country', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('Location', 'eastoutside', 'FontSize', 8);
xlim([0, 600]);
ylim([0, 4.5]);

% Save figure
saveas(gcf, 'inc.png');
hold off;