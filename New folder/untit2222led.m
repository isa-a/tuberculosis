% MATLAB Script to Fill Missing LTBI Values Using Scaling Factor
% This script fills in "Data suppressed" and "[x]" values in true_ltbi.xlsx
% using values from LTBI_prevalence.xlsx multiplied by scaling factor 1.37

clear; clc;

% ==================== CONFIGURATION ====================
scaling_factor = 1.37;  % Median ratio from analysis
true_ltbi_file = 'true ltbi.xlsx';
ltbi_prev_file = 'ghdx.xlsx';
output_file = 'true_ltbi_filled.xlsx';

fprintf('Starting LTBI data filling process...\n\n');

% ==================== LOAD DATA ====================
fprintf('Loading data files...\n');

% Load true_ltbi file
[~, ~, true_ltbi_raw] = xlsread(true_ltbi_file);
true_ltbi_headers = true_ltbi_raw(1, :);  % First row is headers
true_ltbi_data = true_ltbi_raw(2:end, :); % Rest is data

% Load ltbi_prevalence file  
[~, ~, ltbi_prev_raw] = xlsread(ltbi_prev_file);
ltbi_prev_headers = ltbi_prev_raw(1, :);  % First row is headers
ltbi_prev_data = ltbi_prev_raw(2:end, :); % Rest is data

fprintf('  true_ltbi: %d countries loaded\n', size(true_ltbi_data, 1));
fprintf('  ltbi_prevalence: %d countries loaded\n\n', size(ltbi_prev_data, 1));

% ==================== SETUP ====================
% Years to process (2019-2023 are in both files, 2024 only in true_ltbi)
years_both = [2019, 2020, 2021, 2022, 2023];
year_2024_only = 2024;

% Find column indices for years in true_ltbi
year_cols_true = struct();
for i = 1:length(true_ltbi_headers)
    if isnumeric(true_ltbi_headers{i})
        year_cols_true.(sprintf('y%d', true_ltbi_headers{i})) = i;
    end
end

% Find column indices for years in ltbi_prevalence
year_cols_prev = struct();
for i = 1:length(ltbi_prev_headers)
    year_str = ltbi_prev_headers{i};
    if ~ischar(year_str) || strcmp(year_str, 'Country')
        continue;
    end
    year_num = str2double(year_str);
    if ~isnan(year_num)
        year_cols_prev.(sprintf('y%d', year_num)) = i;
    end
end

% Create country mapping (case-insensitive)
fprintf('Creating country mapping...\n');
country_map = containers.Map();
for i = 1:size(ltbi_prev_data, 1)
    country_name = ltbi_prev_data{i, 1};
    if ischar(country_name)
        country_map(lower(strtrim(country_name))) = i;
    end
end
fprintf('  %d countries from prevalence file mapped\n\n', country_map.Count);

% ==================== FILL MISSING VALUES ====================
fprintf('Filling missing values using scaling factor %.2f...\n\n', scaling_factor);

filled_count = 0;
not_filled_count = 0;

% Process each country in true_ltbi
for i = 1:size(true_ltbi_data, 1)
    country_name_true = true_ltbi_data{i, 1};
    
    if ~ischar(country_name_true)
        continue;
    end
    
    country_key = lower(strtrim(country_name_true));
    
    % Check if this country exists in prevalence file
    if ~country_map.isKey(country_key)
        continue;  % Skip if not in prevalence file
    end
    
    prev_row_idx = country_map(country_key);
    
    % Process each year (2019-2023)
    for year = years_both
        year_field = sprintf('y%d', year);
        
        % Get column indices
        true_col = year_cols_true.(year_field);
        prev_col = year_cols_prev.(year_field);
        
        % Get current value from true_ltbi
        current_value = true_ltbi_data{i, true_col};
        
        % Check if value needs to be filled
        needs_filling = false;
        
        if ischar(current_value)
            % Check for "Data suppressed" or "[x]"
            if contains(lower(current_value), 'suppressed') || strcmp(current_value, '[x]')
                needs_filling = true;
            end
        elseif isnan(current_value) || isempty(current_value)
            needs_filling = true;
        end
        
        % Fill if needed
        if needs_filling
            % Get prevalence value
            prev_value = ltbi_prev_data{prev_row_idx, prev_col};
            
            if isnumeric(prev_value) && ~isnan(prev_value)
                % Calculate estimated value
                estimated_value = prev_value * scaling_factor;
                
                % Fill in the value
                true_ltbi_data{i, true_col} = estimated_value;
                
                filled_count = filled_count + 1;
                
                % Print progress every 10 fills
                if mod(filled_count, 10) == 0
                    fprintf('  Filled %d values...\n', filled_count);
                end
            else
                not_filled_count = not_filled_count + 1;
            end
        end
    end
end

fprintf('\n');
fprintf('Processing complete!\n');
fprintf('  Values filled: %d\n', filled_count);
fprintf('  Values not filled (no prevalence data): %d\n', not_filled_count);
fprintf('  Note: 2024 was not filled (no prevalence data for 2024)\n\n');

% ==================== SAVE OUTPUT ====================
fprintf('Saving filled data to: %s\n', output_file);

% Combine headers and data
output_data = [true_ltbi_headers; true_ltbi_data];

% Write to Excel file
xlswrite(output_file, output_data);

fprintf('File saved successfully!\n\n');

% ==================== SUMMARY REPORT ====================
fprintf('==============================================\n');
fprintf('SUMMARY REPORT\n');
fprintf('==============================================\n\n');

fprintf('Input files:\n');
fprintf('  - %s\n', true_ltbi_file);
fprintf('  - %s\n', ltbi_prev_file);
fprintf('\nOutput file:\n');
fprintf('  - %s\n\n', output_file);

fprintf('Scaling factor used: %.2f\n\n', scaling_factor);

fprintf('Formula applied:\n');
fprintf('  Estimated_True_LTBI = LTBI_Prevalence × %.2f\n\n', scaling_factor);

fprintf('Years processed: 2019, 2020, 2021, 2022, 2023\n');
fprintf('Year 2024: NOT filled (no prevalence data available)\n\n');

fprintf('Values filled: %d\n', filled_count);
fprintf('Values unchanged: Existing numeric values preserved\n\n');

fprintf('==============================================\n');
fprintf('Process complete!\n');
fprintf('==============================================\n');