clear; clc; close all;

% Standalone reproduction of Figure 2 from Pig_Result_revision.m.
% The input CSV is a frozen copy of the final 4-stage matrices used by the
% source script, after its subject deletion, zero-to-NaN conversion, INVOS
% baseline fallback, and manual INVOS replacements.

script_dir = fileparts(mfilename('fullpath'));
input_file = fullfile(script_dir, 'figure2_stage_values.csv');

T = readtable(input_file, 'TextType', 'string');
devices = ["NIRSIT-X", "INVOS"];
stage_names = {'Baseline', 'T1', 'T2', 'T3'};

all_stats = table();
fig = figure('Color', 'w', 'Position', [100 100 1500 753], 'Visible', 'off');
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for d = 1:numel(devices)
    Td = T(T.Device == devices(d), :);
    X = Td{:, stage_names};

    stats = run_stagewise_wilcoxon_with_holm(X, devices(d));
    all_stats = [all_stats; stats]; %#ok<AGROW>

    ax = nexttile;
    draw_publication_boxplot(ax, X, stats);
    ylabel(ax, devices(d) + " (%)", 'FontSize', 18, 'FontWeight', 'bold');
end

writetable(all_stats, fullfile(script_dir, 'figure2_stats.csv'));

display_table = table( ...
    all_stats.Device, all_stats.Comparison, ...
    compose('%.3f', all_stats.Raw_p), ...
    compose('%.3f', all_stats.HolmAdjusted_p), ...
    compose('%.3f', all_stats.Z_value), ...
    compose('%.3f', all_stats.EffectSize_r), ...
    'VariableNames', {'Device','Comparison','Raw_p','HolmAdjusted_p','Z_value','EffectSize_r'});
writetable(display_table, fullfile(script_dir, 'figure2_display.csv'));

exportgraphics(fig, fullfile(script_dir, 'figure2_reproduced.png'), 'Resolution', 300);
close(fig);


fprintf('Figure 2 reproduction complete: %s\n', script_dir);

function ResultsTable = run_stagewise_wilcoxon_with_holm(X, device_name)
    comp_names = ["Baseline vs T1"; "Baseline vs T2"; "Baseline vs T3"];
    raw_p = nan(3,1);
    z_vals = nan(3,1);
    r_vals = nan(3,1);
    n_rows = repmat(size(X,1), 3, 1);
    n_valid = nan(3,1);
    n_nonzero = nan(3,1);
    W_vals = nan(3,1);

    baseline = X(:,1);
    for k = 1:3
        y = X(:,k+1);
        valid = ~isnan(baseline) & ~isnan(y);
        x = baseline(valid);
        y = y(valid);
        n_valid(k) = numel(x);

        [p, ~, stats] = signrank(x, y, 'method', 'approximate');
        raw_p(k) = p;
        W_vals(k) = stats.signedrank;

        differences = y - x;
        differences = differences(differences ~= 0);
        n = numel(differences);
        n_nonzero(k) = n;

        if n < 1
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
        repmat(string(device_name),3,1), comp_names, n_rows, n_valid, ...
        n_nonzero, W_vals, z_vals, r_vals, raw_p, adjusted_p, ...
        'VariableNames', {'Device','Comparison','N_rows','N_valid_pairs', ...
        'N_nonzero','SignedRank_W','Z_value','EffectSize_r','Raw_p','HolmAdjusted_p'});
end

function adjusted_p = holm_bonferroni(raw_p)
    adjusted_p = nan(size(raw_p));
    valid = ~isnan(raw_p);
    p_valid = raw_p(valid);
    [p_sorted, order] = sort(p_valid);
    m = numel(p_sorted);
    adjusted_sorted = nan(size(p_sorted));
    for i = 1:m
        adjusted_sorted(i) = min(1, (m-i+1) * p_sorted(i));
    end
    for i = 2:m
        adjusted_sorted(i) = max(adjusted_sorted(i), adjusted_sorted(i-1));
    end
    adjusted_valid = nan(size(p_valid));
    adjusted_valid(order) = adjusted_sorted;
    adjusted_p(valid) = adjusted_valid;
end

function draw_publication_boxplot(ax, X, stats)
    boxplot(ax, X, 'Symbol', '', 'Widths', 0.55);
    hold(ax, 'on');
    ax.Box = 'off';
    ax.LineWidth = 1.8;
    ax.FontSize = 18;
    ax.FontName = 'Arial';
    ax.FontWeight = 'bold';
    ax.TickDir = 'out';
    ax.TickLength = [0.025 0.025];
    ax.XLim = [0.5 4.5];
    ax.YLim = [0 70];
    ax.YTick = 0:10:70;
    ax.XTick = 1:4;
    ax.XTickLabel = {'Baseline', 'T1', 'T2', 'T3'};
    ax.XAxis.FontAngle = 'italic';
    pbaspect(ax, [1 1 1]);

    set(findobj(ax, 'Tag', 'Box'), 'Color', [0 0 0], 'LineWidth', 1.6);
    set(findobj(ax, 'Tag', 'Median'), 'Color', [0.92 0.42 0.42], 'LineWidth', 1.4);
    set(findobj(ax, 'Tag', 'Whisker'), 'Color', [0.55 0.55 0.55], 'LineStyle', '--', 'LineWidth', 1.1);
    set(findobj(ax, 'Tag', 'Upper Whisker'), 'Color', [0.55 0.55 0.55], 'LineStyle', '--', 'LineWidth', 1.1);
    set(findobj(ax, 'Tag', 'Lower Whisker'), 'Color', [0.55 0.55 0.55], 'LineStyle', '--', 'LineWidth', 1.1);

    significant = stats(stats.HolmAdjusted_p < 0.05, :);
    base_y = 54;
    step_y = 4.8;
    for k = 1:height(significant)
        x2 = find(stats.Comparison == significant.Comparison(k), 1) + 1;
        add_sig_bar(ax, 1, x2, base_y + (k-1)*step_y);
    end

    text(ax, 0.94, 0.26, sprintf('N = %d', size(X,1)), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontSize', 22, 'FontAngle', 'italic');
    text(ax, 0.94, 0.16, 'Wilcoxon signed-rank test with Holm correction', ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'FontSize', 18, 'FontAngle', 'italic');
    text(ax, 0.94, 0.08, 'Adjusted p-values: * p < 0.05', ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'FontSize', 18, 'FontAngle', 'italic');
end

function add_sig_bar(ax, x1, x2, y)
    bar_h = 1.5;
    plot(ax, [x1 x1 x2 x2], [y y+bar_h y+bar_h y], ...
        'Color', [0.10 0.20 0.32], 'LineWidth', 1.8);
    text(ax, mean([x1 x2]), y+bar_h+0.25, '*', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 14, 'FontWeight', 'bold');
end
