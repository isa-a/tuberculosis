% Sum "Decisions" by Year (2011–2025) for a chosen country,
% keeping Applicant type = ALL (no filter) and Case outcome = Issued only.
% Reads Book3.xlsx and returns a table with ALL years 2011–2025 (missing -> 0).

%% -------- USER SETTINGS --------
filePath  = "Book3.xlsx";   % change if needed
sheetName = 1;              % or e.g. "Sheet1"
country   = "Angola";  % choose country here
yearMin   = 2011;
yearMax   = 2025;

%% -------- READ DATA (preserve original headers) --------
T = readtable(filePath, ...
              "Sheet", sheetName, ...
              "TextType","string", ...
              "VariableNamingRule","preserve");

%% -------- OPTIONAL: check required columns exist --------
reqCols = ["Year","Nationality","Case outcome","Decisions"];
missing = setdiff(reqCols, string(T.Properties.VariableNames));
if ~isempty(missing)
    error("Missing required column(s): %s", strjoin(missing, ", "));
end

%% -------- CLEAN / COERCE TYPES --------
% Year -> numeric
if ~isnumeric(T.Year)
    T.Year = str2double(string(T.Year));
end

% Decisions -> numeric (handles commas like "3,576")
if ~isnumeric(T.Decisions)
    T.Decisions = str2double( erase(string(T.Decisions), ",") );
else
    % numeric already; do nothing
end

%% -------- FILTER (Applicant type is ALL, so no filter) --------
keep =  (T.Nationality == country) & ...
        (T.("Case outcome") == "Issued") & ...
        (T.Year >= yearMin) & (T.Year <= yearMax);

Tf = T(keep, :);

%% -------- SUM BY YEAR --------
if isempty(Tf)
    % No rows matched; return zeros for all years
    allYears = (yearMin:yearMax)';
    outFull  = table(allYears, zeros(numel(allYears),1), ...
                     'VariableNames', {'Year','TotalDecisions'});
else
    yearTotals = groupsummary(Tf, "Year", "sum", "Decisions");
    yearTotals.Properties.VariableNames(end) = "TotalDecisions";

    % Ensure all years appear (fill missing with 0)
    allYears = (yearMin:yearMax)';
    outFull  = table(allYears, zeros(numel(allYears),1), ...
                     'VariableNames', {'Year','TotalDecisions'});

    [tf2, loc] = ismember(yearTotals.Year, outFull.Year);
    outFull.TotalDecisions(loc(tf2)) = yearTotals.TotalDecisions(tf2);
end

%% -------- OUTPUT --------
disp(outFull);

%% -------- OPTIONAL: write to Excel --------
% outFile = "issued_decisions_by_year_" + replace(country," ","_") + ".xlsx";
% writetable(outFull, outFile, "Sheet", "Totals_2011_2025");
