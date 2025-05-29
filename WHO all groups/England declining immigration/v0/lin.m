years = [2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023];
vals = [797, 772, 752, 825, 886, 682, 773, 1294, 1316];

%y = m*x + b
p = polyfit(years, vals, 1);  %1; slope, 2;intercept

eval = polyval(p, years);


figure;
plot(years, vals, 'ko', 'MarkerFaceColor', 'k'); 
hold on;
plot(years, eval, 'b-', 'LineWidth', 2);         
xlabel('Year');
ylabel('Immigration');
grid on;
