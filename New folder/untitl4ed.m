% Read the CSV
T = readtable('test.csv');

% ---- FORCE numeric and multiply by 10 ----
T.val = double(T.val) * 100;

% ---- Sort countries alphabetically ----
countries = sort(unique(T.location));
years = sort(unique(T.year));

% Preallocate output table
outTable = table(countries, 'VariableNames', {'Country'});

% Fill year columns
for i = 1:numel(years)
    y = years(i);

    % Subset year
    yearRows = T(T.year == y, :);

    % Align rows to alphabetically sorted country list
    [tf, loc] = ismember(countries, yearRows.location);

    col = nan(numel(countries),1);
    col(tf) = yearRows.val(loc(tf));

    outTable.(sprintf('%d', y)) = col;
end

% Write Excel file
writetable(outTable, 'ltbi_prevalence_x10_hi.xlsx');
