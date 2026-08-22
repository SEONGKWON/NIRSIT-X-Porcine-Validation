%% Reproduce Supplementary Table 1 from extracted physiology stage values
% Descriptive statistics use comparison-specific complete pairs.
% The current source code applies Holm correction to four comparisons,
% including T1 versus T2, although the manuscript table displays only the
% three Baseline comparisons.

clear;
clc;

this_dir = fileparts(mfilename('fullpath'));
data_file = fullfile(this_dir, 'physiology_stage_values.csv');
result_file = fullfile(this_dir, 'supple_t1_results.csv');
display_file = fullfile(this_dir, 'supple_t1_display.csv');

opts = detectImportOptions(data_file, 'VariableNamingRule', 'preserve');
opts = setvartype(opts, 'Variable', 'string');
stage_data = readtable(data_file, opts);

required_columns = {'SubjectID','Variable','Baseline','T1','T2','T3'};
assert(all(ismember(required_columns, stage_data.Properties.VariableNames)), ...
    'Input CSV is missing one or more required columns.');

variable_order = ["MBP"; "SBP"; "DBP"; "SpO2"; "HR"; "EtCO2"];
all_results = table();

for v = 1:numel(variable_order)
    variable_name = variable_order(v);
    rows = stage_data.Variable == variable_name;
    assert(sum(rows) == 8, 'Expected eight subject rows for %s.', variable_name);

    X = [stage_data.Baseline(rows), stage_data.T1(rows), ...
         stage_data.T2(rows), stage_data.T3(rows)];
    variable_results = run_comparisons(X, variable_name);
    all_results = [all_results; variable_results]; %#ok<AGROW>
end

writetable(all_results, result_file);

display_rows = all_results.Comparison ~= "T1 vs T2";
shown = all_results(display_rows,:);
display_results = table( ...
    shown.Variable, shown.Comparison, ...
    compose('%.3f [%.3f]', shown.Reference_Median, shown.Reference_IQR), ...
    compose('%.3f [%.3f]', shown.Compared_Median, shown.Compared_IQR), ...
    compose('%.3f', shown.Raw_p), ...
    compose('%.3f', shown.Holm_adjusted_p), ...
    compose('%.3f', shown.Z_value), ...
    compose('%.3f', shown.EffectSize_r), ...
    'VariableNames', {'Variable','Comparison','Baseline_median_IQR', ...
    'Compared_median_IQR','Raw_p','Holm_adjusted_p','Z_value', ...
    'EffectSize_r'});
writetable(display_results, display_file);

disp(' ');
disp('===== Supplementary Table 1: full four-comparison analysis =====');
disp(all_results);
disp(' ');
disp('===== Supplementary Table 1: displayed three comparisons =====');
disp(display_results);


fprintf('Full results: %s\n', result_file);
fprintf('Display table: %s\n', display_file);


function ResultsTable = run_comparisons(X, variable_name)
    % X columns: [Baseline, T1, T2, T3].
    comparison_pairs = [1 2; 1 3; 1 4; 2 3];
    comparison_names = ["Baseline vs T1"; "Baseline vs T2"; ...
                        "Baseline vs T3"; "T1 vs T2"];
    stage_names = ["Baseline"; "T1"; "T2"; "T3"];
    n_comparisons = size(comparison_pairs, 1);

    N_complete = nan(n_comparisons,1);
    N_nonzero = nan(n_comparisons,1);
    ref_median = nan(n_comparisons,1);
    ref_iqr = nan(n_comparisons,1);
    ref_q1 = nan(n_comparisons,1);
    ref_q3 = nan(n_comparisons,1);
    cmp_median = nan(n_comparisons,1);
    cmp_iqr = nan(n_comparisons,1);
    cmp_q1 = nan(n_comparisons,1);
    cmp_q3 = nan(n_comparisons,1);
    raw_p = nan(n_comparisons,1);
    z_value = nan(n_comparisons,1);
    effect_r = nan(n_comparisons,1);

    for c = 1:n_comparisons
        ref_idx = comparison_pairs(c,1);
        cmp_idx = comparison_pairs(c,2);
        x = X(:,ref_idx);
        y = X(:,cmp_idx);

        valid = ~isnan(x) & ~isnan(y);
        x_complete = x(valid);
        y_complete = y(valid);
        N_complete(c) = numel(x_complete);

        qx = prctile(x_complete, [25 75]);
        qy = prctile(y_complete, [25 75]);
        ref_median(c) = median(x_complete);
        ref_q1(c) = qx(1);
        ref_q3(c) = qx(2);
        ref_iqr(c) = qx(2) - qx(1);
        cmp_median(c) = median(y_complete);
        cmp_q1(c) = qy(1);
        cmp_q3(c) = qy(2);
        cmp_iqr(c) = qy(2) - qy(1);

        differences = y_complete - x_complete;
        nonzero = differences ~= 0 & ~isnan(differences);
        x_test = x_complete(nonzero);
        y_test = y_complete(nonzero);
        N_nonzero(c) = numel(x_test);

        if N_nonzero(c) == 0
            raw_p(c) = 1;
            z_value(c) = 0;
            effect_r(c) = 0;
        elseif N_nonzero(c) == 1
            raw_p(c) = NaN;
        else
            [raw_p(c),~,stats] = signrank(x_test, y_test, ...
                'method', 'approximate');
            if isfield(stats, 'zval') && ~isempty(stats.zval)
                z_value(c) = abs(stats.zval);
            elseif raw_p(c) > 0 && raw_p(c) < 1
                z_value(c) = abs(norminv(1 - raw_p(c)/2));
            elseif raw_p(c) == 0
                z_value(c) = Inf;
            else
                z_value(c) = 0;
            end
            effect_r(c) = z_value(c) / sqrt(N_nonzero(c));
        end
    end

    adjusted_p = holm_bonferroni(raw_p);

    ResultsTable = table( ...
        repmat(string(variable_name), n_comparisons, 1), ...
        comparison_names, ...
        stage_names(comparison_pairs(:,1)), ...
        stage_names(comparison_pairs(:,2)), ...
        N_complete, N_nonzero, ...
        ref_median, ref_iqr, ref_q1, ref_q3, ...
        cmp_median, cmp_iqr, cmp_q1, cmp_q3, ...
        raw_p, adjusted_p, z_value, effect_r, ...
        'VariableNames', {'Variable','Comparison','ReferenceStage', ...
        'ComparedStage','N_complete','N_nonzero','Reference_Median', ...
        'Reference_IQR','Reference_Q1','Reference_Q3','Compared_Median', ...
        'Compared_IQR','Compared_Q1','Compared_Q3','Raw_p', ...
        'Holm_adjusted_p','Z_value','EffectSize_r'});
end


function adjusted_p = holm_bonferroni(raw_p)
    adjusted_p = nan(size(raw_p));
    valid = ~isnan(raw_p);
    p = raw_p(valid);
    [p_sorted, sort_idx] = sort(p);
    m = numel(p);
    adjusted_sorted = nan(size(p_sorted));

    for i = 1:m
        adjusted_sorted(i) = (m - i + 1) * p_sorted(i);
    end
    for i = 2:m
        adjusted_sorted(i) = max(adjusted_sorted(i), ...
            adjusted_sorted(i-1));
    end
    adjusted_sorted(adjusted_sorted > 1) = 1;

    p_back = nan(size(p));
    p_back(sort_idx) = adjusted_sorted;
    adjusted_p(valid) = p_back;
end
