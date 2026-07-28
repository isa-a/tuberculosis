% Load data
data = readtable('observed_ratios_only.xlsx');

% Get unique countries
countries = unique(data.Country);
n_countries = length(countries);

% Create figure
figure('Position', [100, 100, 1200, 600]);
hold on;

% Generate unique colors for all 32 countries
% Use multiple colormaps to ensure distinction
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
    
    % Plot with specific color
    plot(country_data.Year, country_data.Ratio, '-o', ...
         'LineWidth', 2, ...
         'MarkerSize', 6, ...
         'Color', colors(i,:), ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:), ...
         'DisplayName', countries{i});
end

% Add horizontal line at 1.37 (median)
yline(1.37, '--r', 'LineWidth', 2, 'DisplayName', 'Median = 1.37');

% Labels and formatting
xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Scaling factor', 'FontSize', 12, 'FontWeight', 'bold');
% title('LTBI Scaling Ratios by Country (2019-2023)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('Location', 'eastoutside', 'FontSize', 8);
xlim([2018.5, 2023.5]);
xticks(2019:2023);  % Only show integer years
xticklabels({'2019', '2020', '2021', '2022', '2023'});  % Force integer labels
ylim([0, 4.5]);

% Save figure
saveas(gcf, 'plot.png');
hold off;