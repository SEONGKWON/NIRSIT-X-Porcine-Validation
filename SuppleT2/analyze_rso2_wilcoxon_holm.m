%% Reproduce Supplementary Table 2 from extracted stage values
% Self-contained extraction of the stage-wise analysis in
% ../Pig_Result_revision.m. Requires MATLAB and Statistics and Machine
% Learning Toolbox for signrank.

clear;
clc;

this_dir = fileparts(mfilename('fullpath'));
data_file = fullfile(this_dir, 'rso2_stage_values.csv');
result_file = fullfile(this_dir, 'supple_t2_results.csv');
display_file = fullfile(this_dir, 'supple_t2_display.csv');

opts = detectImportOptions(data_file, 'VariableNamingRule', 'preserve');
opts = setvartype(opts, 'Device', 'string');
stage_data = readtable(data_file, opts);

required_columns = {'SubjectID','Device','Baseline','T1','T2','T3'};
assert(all(ismember(required_columns, stage_data.Properties.VariableNames)), ...
    'Input CSV is missing one or more required columns.');

device_order = ["NIRSIT-X"; "INVOS"];
all_results = table();

for i = 1:numel(device_order)
    device_name = device_order(i);
    rows = stage_data.Device == device_name;
    assert(any(rows), 'No rows found for device %s.', device_name);

    X = [stage_data.Baseline(rows), stage_data.T1(rows), ...
         stage_data.T2(rows), stage_data.T3(rows)];
    device_results = run_stagewise_wilcoxon_with_holm(X, device_name);
    all_results = [all_results; device_results]; %#ok<AGROW>
end

writetable(all_results, result_file);

display_results = table( ...
    all_results.Device, ...
    all_results.Comparison, ...
    compose('%.3f', all_results.Raw_p), ...
    compose('%.3f', all_results.HolmBonferroni_p), ...
    compose('%.3f', all_results.Z_value), ...
    compose('%.3f', all_results.EffectSize_r), ...
    'VariableNames', { ...
        'Device','Comparison','Raw_p','Holm_adjusted_p', ...
        'Z_value','Effect_size_r'});
writetable(display_results, display_file);

disp(' ');
disp('===== Supplementary Table 2: full analysis precision =====');
disp(all_results);
disp(' ');
disp('===== Supplementary Table 2: manuscript display precision =====');
disp(display_results);


fprintf('Full results: %s\n', result_file);
fprintf('Display table: %s\n', display_file);


function ResultsTable = run_stagewise_wilcoxon_with_holm(X, device_name)
% X is n_subjects x 4: [Baseline, T1, T2, T3].
% This function intentionally reproduces Pig_Result_revision.m:
%   raw p: MATLAB signrank(..., 'method', 'approximate')
%   Z: manual normal approximation with continuity correction
%   r: abs(Z) / sqrt(N_nonzero)

    comp_names = {'Baseline vs T1'; 'Baseline vs T2'; 'Baseline vs T3'};
    raw_p = nan(3,1);
    z_vals = nan(3,1);
    r_vals = nan(3,1);
    n_pairs = nan(3,1);
    W_vals = nan(3,1);

    baseline = X(:,1);

    for k = 1:3
        y = X(:,k+1);

        valid = ~isnan(baseline) & ~isnan(y);
        x_valid = baseline(valid);
        y_valid = y(valid);

        [p, ~, stats] = signrank(x_valid, y_valid, ...
            'method', 'approximate');
        raw_p(k) = p;
        W_vals(k) = stats.signedrank;

        differences = y_valid - x_valid;
        differences_nonzero = differences(differences ~= 0);
        n = numel(differences_nonzero);
        n_pairs(k) = n;

        if n < 1
            z_vals(k) = NaN;
            r_vals(k) = NaN;
            continue;
        end

        mu_W = n * (n + 1) / 4;
        sigma_W = sqrt(n * (n + 1) * (2*n + 1) / 24);

        if stats.signedrank > mu_W
            z = (stats.signedrank - mu_W - 0.5) / sigma_W;
        elseif stats.signedrank < mu_W
            z = (stats.signedrank - mu_W + 0.5) / sigma_W;
        else
            z = 0;
        end

        z_vals(k) = z;
        r_vals(k) = abs(z) / sqrt(n);
    end

    adjusted_p = holm_bonferroni(raw_p);

    ResultsTable = table( ...
        repmat(string(device_name), 3, 1), ...
        string(comp_names), ...
        n_pairs, ...
        W_vals, ...
        z_vals, ...
        r_vals, ...
        raw_p, ...
        adjusted_p, ...
        'VariableNames', { ...
            'Device','Comparison','N_nonzero','SignedRank_W','Z_value', ...
            'EffectSize_r','Raw_p','HolmBonferroni_p'});
end


function adjusted_p = holm_bonferroni(raw_p)
    [p_sorted, sort_idx] = sort(raw_p);
    m = numel(raw_p);
    adjusted_sorted = nan(size(p_sorted));

    for i = 1:m
        adjusted_sorted(i) = min((m - i + 1) * p_sorted(i), 1);
    end

    for i = 2:m
        adjusted_sorted(i) = max(adjusted_sorted(i), ...
            adjusted_sorted(i-1));
    end

    adjusted_p = nan(size(raw_p));
    adjusted_p(sort_idx) = adjusted_sorted;
end
