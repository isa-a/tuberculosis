% Data
% relapse_risks = {'3.2, 14', '1.6, 7', '0.8, 3.5', '0.4, 1.75', '0.2, 0.875'};  
relapse_risks = {'Existing relapse rates', '50% reduction', '75% reduction', '87.5% reduction', '95% reduction'};  


% Incidences for different TPT eff levels
incidences_40 = [2.02, 1.71, 1.57, 1.5, 1];
incidences_60 = [1.24, 0.96, 0.84, 0.78, 0.75];
incidences_80 = [0.45, 0.2, 0.1, 0.05, 0.02];

% Create the plot with three curves
ff=figure('Position', [577,   190 ,   1029 ,732]); lw = 1.5; fs = 14;
% Plot for TPT Eff 0.8 (blue)
plot(incidences_40, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'blue', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue', 'DisplayName', 'TPT Eff 0.8');
hold on;
% Plot for TPT Eff 0.9 (green)
% plot(incidences_60, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'green', 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green', 'DisplayName', 'TPT Eff 0.9');
% % Plot for TPT Eff 1 (red)
% plot(incidences_80, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'DisplayName', 'TPT Eff 1');
% hold off;

% Set x-axis labels as the relapse risks
set(gca, 'XTick', 1:length(relapse_risks), 'XTickLabel', relapse_risks);

% Add labels and title
xlabel('Relapse rates', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Incidence in 2040 (per 100,000)', 'FontWeight', 'bold', 'FontSize', 12);
%title('Incidence vs Relapse Risk at Different TPT Eff Levels', 'FontWeight', 'bold', 'FontSize', 14);

% Set y-axis to start at 0
ylim([0, max(incidences_40) + 0.5]);

% Add a dashed horizontal line at y = 0
yline(0.1, '--', 'LineWidth', 1.5, 'Color', 'k','HandleVisibility', 'off');
yline(1, '--', 'LineWidth', 1.5, 'Color', 'k','HandleVisibility', 'off');

% Add legend
legend({'80% effective TPT', '90% effective TPT', '100% effective TPT'}, 'FontWeight', 'bold', 'FontSize', 12);

% Display grid
grid on;

% Adjust layout
set(gca, 'FontSize', 12, 'FontWeight', 'bold');



%---------------------------------------------------------------vaccine


% Data
% relapse_risks = {'3.2, 14', '1.6, 7', '0.8, 3.5', '0.4, 1.75', '0.2, 0.875'};  
relapse_risks = {'Existing relapse rates', '50% reduction', '75% reduction', '87.5% reduction', '93.75% reduction'};  


% Incidences for different TPT eff levels
incidences_40 = [2.45,2.13,1.98,1.91, 1.87];
incidences_60 = [1.89,1.59,1.45,1.39,1.35];
incidences_80 = [1.36,1.05,0.93,0.87,0.83];

% Create the plot with three curves
ff=figure('Position', [577,   190 ,   1029 ,732]); lw = 1.5; fs = 14;
% Plot for TPT Eff 0.8 (blue)
plot(incidences_40, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'blue', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue', 'DisplayName', 'TPT Eff 0.8');
hold on;
% Plot for TPT Eff 0.9 (green)
plot(incidences_60, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'green', 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green', 'DisplayName', 'TPT Eff 0.9');
% Plot for TPT Eff 1 (red)
plot(incidences_80, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'DisplayName', 'TPT Eff 1');
hold off;

% Set x-axis labels as the relapse risks
set(gca, 'XTick', 1:length(relapse_risks), 'XTickLabel', relapse_risks);

% Add labels and title
xlabel('Relapse rates', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Incidence in 2040 (per 100,000)', 'FontWeight', 'bold', 'FontSize', 12);
%title('Incidence at varying relapse rates and vaccine efficacy, TPT 60% effective', 'FontWeight', 'bold', 'FontSize', 14);

% Set y-axis to start at 0
ylim([0, max(incidences_40) + 0.5]);

% Add a dashed horizontal line at y = 0
yline(0.1, '--', 'LineWidth', 1.5, 'Color', 'k','HandleVisibility', 'off');
yline(1, '--', 'LineWidth', 1.5, 'Color', 'k','HandleVisibility', 'off');

% Add legend
legend({'40% effective vaccine', '60% effective vaccine', '80% effective vaccine'}, 'FontWeight', 'bold', 'FontSize', 12);

% Display grid
grid on;

% Adjust layout
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

%---------------------------------------------------------------vaccine2


% Data
% relapse_risks = {'3.2, 14', '1.6, 7', '0.8, 3.5', '0.4, 1.75', '0.2, 0.875'};  
relapse_risks = {'Existing relapse rates', '50% reduction', '75% reduction', '87.5% reduction', '93.75% reduction'};  


% Incidences for different TPT eff levels
incidences_40 = [0.45,0.21,0.1,0.05,0.024];
incidences_60 = [1.89,1.59,1.45,1.39,1.35];
incidences_80 = [1.36,1.05,0.93,0.87,0.83];

% Create the plot with three curves
ff=figure('Position', [577,   190 ,   1029 ,732]); lw = 1.5; fs = 14;
% Plot for TPT Eff 0.8 (blue)
plot(incidences_40, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'blue', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue', 'DisplayName', 'TPT Eff 0.8');
hold on;
% Plot for TPT Eff 0.9 (green)
% plot(incidences_60, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'green', 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green', 'DisplayName', 'TPT Eff 0.9');
% % Plot for TPT Eff 1 (red)
% plot(incidences_80, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'DisplayName', 'TPT Eff 1');
% hold off;

% Set x-axis labels as the relapse risks
set(gca, 'XTick', 1:length(relapse_risks), 'XTickLabel', relapse_risks);

% Add labels and title
xlabel('Relapse rates', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Incidence in 2040 (per 100,000)', 'FontWeight', 'bold', 'FontSize', 12);
%title('Incidence at varying relapse rates and vaccine efficacy, TPT 60% effective', 'FontWeight', 'bold', 'FontSize', 14);

% Set y-axis to start at 0
ylim([0, max(incidences_40) + 0.5]);

% Add a dashed horizontal line at y = 0
yline(0.1, '--', 'LineWidth', 1.5, 'Color', 'k','HandleVisibility', 'off');
yline(1, '--', 'LineWidth', 1.5, 'Color', 'k','HandleVisibility', 'off');

% Add legend
legend({'40% effective vaccine', '60% effective vaccine', '80% effective vaccine'}, 'FontWeight', 'bold', 'FontSize', 12);

% Display grid
grid on;

% Adjust layout
set(gca, 'FontSize', 12, 'FontWeight', 'bold');