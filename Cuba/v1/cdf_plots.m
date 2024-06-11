load analysis_outputs.mat

figure; histogram(ch_inc)
ax = gca;
ax.FontSize = 14; 
ax.FontWeight = 'bold';
xlabel('Incidence in children (per 100,000 children-years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Frequency', 'FontSize', 14, 'FontWeight', 'bold');


[yy, xx] = ecdf(ch_inc);
%[yy2, xx2] = ecdf(ch_inc2);


figure; hold on;
plot(xx, 1-yy, 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); 
xlim([min(xx), 0.2]);
ylim([0, 1]);
ax = gca;
ax.FontSize = 14; 
ax.FontWeight = 'bold';
xlabel('Incidence in children (per 100,000 children-years)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');




boundary_66 = xx(find(yy >= 0.66, 1));
boundary_33 = xx(find(yy >= 0.33, 1));


line([boundary_66, boundary_66], [0, 1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
line([boundary_33, boundary_33], [0, 1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');


fill([min(xx); xx(xx <= boundary_33); boundary_33], [0; 1-yy(xx <= boundary_33); 0], [1, 0.84, 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Gold'); % Gold
fill([boundary_33; xx(xx > boundary_33 & xx <= boundary_66); boundary_66], [0; 1-yy(xx > boundary_33 & xx <= boundary_66); 0], [0.75, 0.75, 0.75], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Silver'); % Silver
fill([boundary_66; xx(xx > boundary_66); max(xx)], [0; 1-yy(xx > boundary_66); 0], [0.8, 0.5, 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Bronze'); % Bronze
% Add horizontal lines at y = 0.66 and y = 0.33
line([0, boundary_66], [1-0.66, 1-0.66], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
line([0, boundary_33], [1-0.33, 1-0.33], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
% Add filled circle dots where the lines meet
plot(boundary_66, 1-0.66, 'ko', 'MarkerFaceColor', 'k');
plot(boundary_33, 1-0.33, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5); 


xlim([min(xx), 0.18]);
title('Children incidence CDF', 'FontWeight', 'bold');
xlabel('Incidence in children (per 100,000 children-years)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');
legend('Location', 'best', 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold');
hold off;





%%%%%%

% get data 
[yy, xx] = ecdf(ch_inc2);
%[yy2, xx2] = ecdf(ch_inc2);


figure; hold on;
plot(xx, 1-yy, 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); 

% get bounds
boundary_66 = xx(find(yy >= 0.66, 1));
boundary_33 = xx(find(yy >= 0.33, 1));

% plot verticals
line([boundary_66, boundary_66], [0, 1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
line([boundary_33, boundary_33], [0, 1], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

% shade
fill([min(xx); xx(xx <= boundary_33); boundary_33], [0; 1-yy(xx <= boundary_33); 0], [1, 0.84, 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Gold'); % Gold
fill([boundary_33; xx(xx > boundary_33 & xx <= boundary_66); boundary_66], [0; 1-yy(xx > boundary_33 & xx <= boundary_66); 0], [0.75, 0.75, 0.75], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Silver'); % Silver
fill([boundary_66; xx(xx > boundary_66); max(xx)], [0; 1-yy(xx > boundary_66); 0], [0.8, 0.5, 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Bronze'); % Bronze

%  x-axis limit
%xlim([min(xx), 0.18]);

title('Children incidence CDF', 'FontWeight', 'bold');
xlabel('Incidence in children (cases per 100,000 children)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');

legend('Location', 'best', 'FontWeight', 'bold');

set(gca, 'FontWeight', 'bold');

hold off;



[yy1, xx1] = ecdf(ch_inc);
[yy2, xx2] = ecdf(ch_inc2);


figure; hold on;
plot(xx1, 1-yy1, 'Color', [65/255, 80/255, 225/255], 'LineWidth', 2); %  blue
plot(xx2, 1-yy2, 'Color', [255/255, 0, 0], 'LineWidth', 2); %  red


ax = gca;
ax.FontSize = 14; 
ax.FontWeight = 'bold';
xlabel('Incidence in children (cases per 100,000 children)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontSize', 14, 'FontWeight', 'bold');
title('Thresholds for reaching elimination vs pre-elimination', 'FontSize', 16, 'FontWeight', 'bold');

hold off;


%tornado plot


vec = pcm(end, 1:end-1);


param_names = {'Contact rate', 'Secular trend in contact rate', 'Per-capita rate of treatment initiation', 'Child mortality', 'Relapse risk', 'Ageing', 'Relative infectivity, children vs adults', 'Intergenerational mixing', 'Baseline incidence'};


[~, sorted_indices] = sort(abs(vec), 'descend');
sorted_vec = vec(sorted_indices);


sorted_param_names = param_names(sorted_indices);

sorted_vec = flipud(sorted_vec);
sorted_param_names = flipud(sorted_param_names);


figure;
barh(sorted_vec, 'red');


set(gca, 'YDir', 'reverse');


xlabel('Partial Correlation', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Parameters', 'FontSize', 14, 'FontWeight', 'bold');
title('Partial Correlation Ordered by Absolute Value', 'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'YTick', 1:length(sorted_param_names), 'YTickLabel', sorted_param_names);


set(gca, 'FontWeight', 'bold');




%%%%%%

% re do no shading with bounds




[yy, xx] = ecdf(ch_inc);
%[yy2, xx2] = ecdf(ch_inc2);

figure; hold on;
plot(xx, 1-yy, 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); 
xlim([min(xx), 0.2]);
ylim([0, 1]);
ax = gca;
ax.FontSize = 14; 
ax.FontWeight = 'bold';
title('Children incidence CDF', 'FontWeight', 'bold');
xlabel('Incidence in children (cases per 100,000 children)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');

boundary_66 = xx(find(yy >= 0.66, 1));
boundary_33 = xx(find(yy >= 0.33, 1));

line([boundary_66, boundary_66], [0, 1-0.66], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

% Add horizontal lines at y = 0.66 and y = 0.33
line([0, boundary_66], [1-0.66, 1-0.66], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');

% Add filled circle dots where the lines meet
plot(boundary_66, 1-0.66, 'ko', 'MarkerFaceColor', 'k');

xlim([min(xx), 0.18]);
title('Children incidence CDF', 'FontWeight', 'bold');
xlabel('Incidence in children (per 100,000 children-years)', 'FontWeight', 'bold');
ylabel('Probability that country has reached elimination', 'FontWeight', 'bold');
legend('Location', 'best', 'FontWeight', 'bold');
set(gca, 'FontWeight', 'bold');
hold off;
