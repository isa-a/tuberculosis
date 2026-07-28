


filePath  = "Book3.xlsx";  
sheetName = 1;
yearMin   = 2011;
yearMax   = 2025;

countries = [

"Albania"
"Algeria"
"Armenia"
"Azerbaijan"
"Belarus"
"Benin"
"Bolivia, Plurinational State of"
"Bosnia and Herzegovina"
"Burkina Faso"
"Burundi"
"Cabo Verde"
"Chad"
"China"
"Colombia"
"Comoros"
"Cuba"
"Dominican Republic"
"Ecuador"
"Egypt"
"El Salvador"
"Eritrea"
"Ethiopia"
"Fiji"
"Georgia"
"Ghana"
"Guatemala"
"Honduras"
"Iran, Islamic Republic of"
"Iraq"
"Jamaica"
"Jordan"
"Kazakhstan"
"Kyrgyzstan"
"Lao, People's Democratic Republic"
"Lebanon"
"Libya"
"Malawi"
"Mali"
"Mauritania"
"Moldova, Republic of"
"Montenegro"
"Morocco"
"Niger"
"North Macedonia"
"Occupied Palestinian territory, including east Jerusalem"
"Rwanda"
"Senegal"
"Serbia"
"Sri Lanka"
"Sudan"
"Syrian Arab Republic"
"Tajikistan"
"Togo"
"Tunisia"
"Turkey"
"Turkmenistan"
"Ukraine"
"Uzbekistan"
"Vanuatu"
"Venezuela, Bolivarian Republic of"
"Yemen"

];



appTypes = ["All"];

%read
T = readtable(filePath, ...
              "Sheet", sheetName, ...
              "TextType","string", ...
              "VariableNamingRule","preserve");

%check cols
reqCols = ["Year","Nationality","Case outcome","Applicant type","Decisions"];
missing = setdiff(reqCols, string(T.Properties.VariableNames));
if ~isempty(missing)
    error("Missing required column(s): %s", strjoin(missing, ", "));
end


T.Year = str2double(string(T.Year));
T.Decisions = str2double( erase(string(T.Decisions), ",") );

%excel struct
years = (yearMin:yearMax)';
nY = numel(years);
nC = numel(countries);
resultMat = zeros(nC, nY);


for i = 1:nC

    country = countries(i);

    keep =  (T.Nationality == country) & ...
            (T.("Case outcome") == "Issued") & ...
            ismember(T.("Applicant type"), appTypes) & ...
            (T.Year >= yearMin) & (T.Year <= yearMax);

    Tf = T(keep, :);

    if ~isempty(Tf)
        tmp = groupsummary(Tf, "Year", "sum", "Decisions");
        tmp.Properties.VariableNames(end) = "TotalDecisions";

        [tf2, loc] = ismember(tmp.Year, years);
        resultMat(i, loc(tf2)) = tmp.TotalDecisions(tf2);
    end
end

%write
outTable = array2table(resultMat, 'VariableNames', "Y" + string(years));
outTable = addvars(outTable, countries, 'Before', 1, 'NewVariableNames', "Country");


disp(outTable)
writetable(outTable, "visastbl_lo.xlsx");
