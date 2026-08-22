clear; clc; close all;

% Standalone reproduction of the LMM panels in LMM_AI_IN.m (Figure 3).
% Column orientation intentionally follows the source calculation:
% Y = workbook column 2 (actual rSO2), z = workbook column 3 (actual SjvO2).

script_dir = fileparts(mfilename('fullpath'));
input_file = fullfile(script_dir, 'figure_input_data.csv');

Tinput = readtable(input_file, 'TextType', 'string');
Tinput = Tinput(Tinput.IncludeInAnalysis == 1, :);
devices = ["NIRSIT-X", "INVOS"];

results = cell(numel(devices),1);
result_table = table();

for d = 1:numel(devices)
    Td = Tinput(Tinput.Device == devices(d), :);
    [results{d}, row] = run_lmm_analysis(Td, devices(d));
    result_table = [result_table; row]; %#ok<AGROW>
end

writetable(result_table, fullfile(script_dir, 'figure3_results.csv'));

display_table = table( ...
    result_table.Device, result_table.NObservations, result_table.NSubjects, ...
    compose('%.3f', result_table.Intercept), ...
    compose('%.3f', result_table.Slope), ...
    compose('%.4g', result_table.SlopeP), ...
    'VariableNames', {'Device','NObservations','NSubjects','Intercept','Slope', ...
    'SlopeP'});
writetable(display_table, fullfile(script_dir, 'figure3_display.csv'));

fig = figure('Color', 'w', 'Position', [100 100 1400 650], 'Visible', 'off');
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for d = 1:numel(devices)
    ax = nexttile;
    draw_lmm_panel(ax, results{d});
end
exportgraphics(fig, fullfile(script_dir, 'figure3_reproduced.png'), 'Resolution', 300);
close(fig);


fprintf('Figure 3 reproduction complete: %s\n', script_dir);

function [result, output_row] = run_lmm_analysis(Td, device)
    Y = Td.rSO2;       % source Excel column 2
    z = Td.SjvO2;      % source Excel column 3
    subject_numeric = Td.SubjectID;
    Subject = categorical(subject_numeric);
    T = table(Y(:), z(:), Subject(:), 'VariableNames', {'Y','z','Subject'});

    lme = fitlme(T, 'Y ~ z + (1|Subject)');
    coefficients = lme.Coefficients;
    beta = fixedEffects(lme);
    intercept = beta(1);
    slope = beta(2);

    x_fit = linspace(min(z), max(z), 100)';
    fake_subject = categorical(repmat(subject_numeric(1), size(x_fit)));
    T_predict = table(x_fit, fake_subject, 'VariableNames', {'z','Subject'});
    [y_fit, y_ci] = predict(lme, T_predict);

    random_parameters = covarianceParameters(lme);
    tau_00 = random_parameters{1,1};
    sigma2 = lme.MSE;
    ICC = tau_00/(tau_00+sigma2);

    subject_levels = categories(T.Subject);
    n_subjects = numel(subject_levels);

    output_row = table( ...
        string(device), lme.NumObservations, n_subjects, ...
        coefficients.Estimate(1), coefficients.SE(1), coefficients.pValue(1), ...
        coefficients.Estimate(2), coefficients.SE(2), coefficients.pValue(2), ...
        sigma2, tau_00, ICC, ...
        'VariableNames', {'Device','NObservations','NSubjects', ...
        'Intercept','InterceptSE','InterceptP','Slope','SlopeSE','SlopeP', ...
        'ResidualVariance','RandomInterceptVariance','ICC'});

    result = struct('device',device,'Y',Y,'z',z,'subject',subject_numeric, ...
        'x_fit',x_fit,'y_fit',y_fit,'y_ci',y_ci,'intercept',intercept, ...
        'slope',slope,'slope_p',coefficients.pValue(2));
end

function draw_lmm_panel(ax, result)
    hold(ax, 'on');
    subjects = unique(result.subject);
    colors = lines(numel(subjects));
    for i = 1:numel(subjects)
        idx = result.subject == subjects(i);
        scatter(ax, result.z(idx), result.Y(idx), 50, 'filled', ...
            'MarkerFaceColor', colors(i,:), 'HandleVisibility', 'off');
    end
    fill(ax, [result.x_fit; flipud(result.x_fit)], ...
        [result.y_ci(:,1); flipud(result.y_ci(:,2))], 'r', ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, result.x_fit, result.y_fit, 'r-', 'LineWidth', 2, ...
        'HandleVisibility', 'off');
    xlabel(ax, 'rSO2 (%)');
    ylabel(ax, 'SjvO2 (%)');
    title(ax, sprintf('%s LMM Fit: SjvO2 = %.3f * rSO2 + %.3f (p = %.3g)', ...
        result.device, result.slope, result.intercept, result.slope_p));
    grid(ax, 'on');
end
