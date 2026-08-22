%% Reproduce manuscript Table 1 from extracted physiology stage values
% Uses all available observations separately at each stage.

clear;
clc;

this_dir = fileparts(mfilename('fullpath'));
data_file = fullfile(this_dir, 'physiology_stage_values.csv');
result_file = fullfile(this_dir, 'table1_results.csv');
display_file = fullfile(this_dir, 'table1_display.csv');

opts = detectImportOptions(data_file, 'VariableNamingRule', 'preserve');
opts = setvartype(opts, 'Variable', 'string');
stage_data = readtable(data_file, opts);

required_columns = {'SubjectID','Variable','Baseline','T1','T2','T3'};
assert(all(ismember(required_columns, stage_data.Properties.VariableNames)), ...
    'Input CSV is missing one or more required columns.');

variable_order = ["MBP"; "SBP"; "DBP"; "SpO2"; "HR"; "EtCO2"];
stage_names = ["Baseline"; "T1"; "T2"; "T3"];
all_results = table();

for v = 1:numel(variable_order)
    variable_name = variable_order(v);
    rows = stage_data.Variable == variable_name;
    assert(sum(rows) == 8, 'Expected eight subject rows for %s.', variable_name);

    X = [stage_data.Baseline(rows), stage_data.T1(rows), ...
         stage_data.T2(rows), stage_data.T3(rows)];

    for s = 1:numel(stage_names)
        values = X(:,s);
        values = values(~isnan(values));
        assert(~isempty(values), 'No observations for %s at %s.', ...
            variable_name, stage_names(s));

        q = prctile(values, [25 75]);
        result_row = table( ...
            variable_name, stage_names(s), numel(values), ...
            mean(values), std(values, 0), median(values), ...
            q(2) - q(1), q(1), q(2), ...
            'VariableNames', {'Variable','Stage','N','Mean','SD', ...
            'Median','IQR','Q1','Q3'});
        all_results = [all_results; result_row]; %#ok<AGROW>
    end
end

writetable(all_results, result_file);

display_results = table(variable_order, ...
    'VariableNames', {'Variable'});
for s = 1:numel(stage_names)
    rows = all_results.Stage == stage_names(s);
    display_results.(stage_names(s)) = compose('%.2f [%.2f]', ...
        all_results.Median(rows), all_results.IQR(rows));
    display_results.(stage_names(s) + "_N") = all_results.N(rows);
end
writetable(display_results, display_file);

disp(' ');
disp('===== Manuscript Table 1: full analysis precision =====');
disp(all_results);
disp(' ');
disp('===== Manuscript Table 1: two-decimal display =====');
disp(display_results);


fprintf('Full results: %s\n', result_file);
fprintf('Display table: %s\n', display_file);
