filePath  = "Book3.xlsx";
sheetName = 1;
yearMin   = 2011;
yearMax   = 2025;
countries = [
"Afghanistan"
"Angola"
"Bangladesh"
"Benin"
"Bhutan"
"Botswana"
"Burkina Faso"
"Burundi"
"Cape Verde"
"Cambodia"
"Cameroon"
"Central African Republic"
"Chad"
"Comoros"
"Congo"
"Congo (Democratic Republic)"
"Djibouti"
"Equatorial Guinea"
"Eritrea"
"Eswatini"
"Ethiopia"
"Gabon"
"Gambia, The"
"Ghana"
"Greenland"
"Guinea"
"Guinea-Bissau"
"Haiti"
"India"
"Indonesia"
"Ivory Coast"
"Kenya"
"Kiribati"
"Lesotho"
"Liberia"
"Madagascar"
"Malawi"
"Mali"
"Marshall Islands"
"Mauritania"
"Mauritius"
"Mongolia"
"Mozambique"
"Myanmar (Burma)"
"Namibia"
"Nauru"
"Nepal"
"Niger"
"Nigeria"
"North Korea"
"Pakistan"
"Papua New Guinea"
"Peru"
"Philippines"
"Rwanda"
"Sao Tome and Principe"
"Senegal"
"Seychelles"
"Sierra Leone"
"Somalia"
"South Africa"
"South Sudan"
"Sudan"
"Tanzania"
"Thailand"
"East Timor"
"Togo"
"Uganda"
"Vietnam"
"Zambia"
"Zimbabwe"
];
appTypes = ["Main Applicant", "Dependant"];

% read
T = readtable(filePath, ...
    "Sheet", sheetName, ...
    "TextType","string", ...
    "VariableNamingRule","preserve");

% check cols
reqCols = ["Year","Nationality","Applicant type","Case outcome","Decisions"];
missing = setdiff(reqCols, string(T.Properties.VariableNames));
if ~isempty(missing)
    error("Missing required column(s): %s", strjoin(missing, ", "));
end

T.Year = str2double(string(T.Year));
T.Decisions = str2double(erase(string(T.Decisions), ","));

% build result matrix
years = (yearMin:yearMax)';
nY = numel(years);
nC = numel(countries);
resultMat = zeros(nC, nY);

for i = 1:nC
    country = countries(i);
    keep =  (T.Nationality == country) & ...
        (T.("Case outcome") == "Issued") & ...
        ismember(T.("Applicant type"), appTypes) & ...
        (T.Year >= yearMin) & (T.Year <= yearMax) & ...
        (T.("Visa type group") ~= "Visitor") & ...
        (T.("Visa type subgroup") ~= "Transit");
    Tf = T(keep, :);
    if ~isempty(Tf)
        tmp = groupsummary(Tf, "Year", "sum", "Decisions");
        tmp.Properties.VariableNames(end) = "TotalDecisions";
        [tf2, loc] = ismember(tmp.Year, years);
        resultMat(i, loc(tf2)) = tmp.TotalDecisions(tf2);
    end
end

% write output
outTable = array2table(resultMat, 'VariableNames', "Y" + string(years));
outTable = addvars(outTable, countries, 'Before', 1, 'NewVariableNames', "Country");
disp(outTable)
writetable(outTable, "hitbvisas.xlsx");