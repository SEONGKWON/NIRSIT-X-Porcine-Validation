clear; clc; close all;

% Standalone reproduction of Four_axis_relative_3.m (Figure 4).
% It preserves the source script's subject-specific pseudo-baseline rule and
% its OR exclusion rule: abs(x)>threshold OR abs(y)>threshold.

script_dir = fileparts(mfilename('fullpath'));
input_file = fullfile(script_dir, 'figure_input_data.csv');

T = readtable(input_file, 'TextType', 'string');
T = T(T.IncludeInAnalysis == 1, :);
devices = ["NIRSIT-X", "INVOS"];
thresholds = [5 10];

all_points = table();
all_results = table();
fig = figure('Color', 'w', 'Position', [100 100 1400 650], 'Visible', 'off');
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for d = 1:numel(devices)
    Td = T(T.Device == devices(d), :);
    [r_delta, s_delta, is_baseline] = deltas_from_pseudo_baseline( ...
        Td.rSO2, Td.SjvO2, Td.SubjectID);

    point_table = table( ...
        Td.Device, Td.SourceRow, Td.SubjectID, Td.rSO2, Td.SjvO2, ...
        is_baseline, r_delta, s_delta, ...
        'VariableNames', {'Device','SourceRow','SubjectID','rSO2','SjvO2', ...
        'IsPseudoBaseline','rSO2RelativeChange','SjvO2RelativeChange'});
    all_points = [all_points; point_table]; %#ok<AGROW>

    results = calculate_concordance(devices(d), r_delta, s_delta, ...
        is_baseline, thresholds);
    all_results = [all_results; results]; %#ok<AGROW>

    ax = nexttile;
    draw_four_quadrant_panel(ax, r_delta, s_delta, is_baseline, ...
        results, devices(d));
end

writetable(all_points, fullfile(script_dir, 'figure4_point_data.csv'));
writetable(all_results, fullfile(script_dir, 'figure4_results.csv'));

display_table = table( ...
    all_results.Device, all_results.ExclusionThreshold, ...
    all_results.ConcordantCount, all_results.NIncluded, ...
    compose('%.1f%%', 100*all_results.Concordance), ...
    compose('%.1f%%', 100*all_results.CI95_Lower), ...
    compose('%.1f%%', 100*all_results.CI95_Upper), ...
    'VariableNames', {'Device','ExclusionThreshold','ConcordantCount', ...
    'NIncluded','Concordance','CI95_Lower','CI95_Upper'});
writetable(display_table, fullfile(script_dir, 'figure4_display.csv'));

exportgraphics(fig, fullfile(script_dir, 'figure4_reproduced.png'), 'Resolution', 300);
close(fig);


fprintf('Figure 4 reproduction complete: %s\n', script_dir);

function [dx, dy, is_baseline] = deltas_from_pseudo_baseline(r, s, subject)
    subjects = unique(subject, 'stable');
    dx = nan(size(r));
    dy = nan(size(s));
    is_baseline = false(size(r));

    for i = 1:numel(subjects)
        idx = find(subject == subjects(i));
        rg = r(idx);
        sg = s(idx);
        median_r = median(rg, 'omitnan');
        median_s = median(sg, 'omitnan');
        distance = hypot(rg-median_r, sg-median_s);
        distance(~isfinite(distance)) = Inf;
        if all(~isfinite(distance))
            local_baseline = 1;
        else
            [~, local_baseline] = min(distance);
        end
        baseline_idx = idx(local_baseline);
        rb = r(baseline_idx);
        sb = s(baseline_idx);

        if isfinite(rb) && abs(rb) > eps
            dx(idx) = 100*(r(idx)-rb)./rb;
            dx(baseline_idx) = 0;
        end
        if isfinite(sb) && abs(sb) > eps
            dy(idx) = 100*(s(idx)-sb)./sb;
            dy(baseline_idx) = 0;
        end
        is_baseline(baseline_idx) = true;
    end
end

function results = calculate_concordance(device, dx, dy, is_baseline, thresholds)
    valid = ~isnan(dx) & ~isnan(dy);
    x = dy(valid);  % source x-axis: SjvO2 relative change
    y = dx(valid);  % source y-axis: rSO2 relative change
    baseline = is_baseline(valid);

    n_thresholds = numel(thresholds);
    concordant_count = nan(n_thresholds,1);
    n_included = nan(n_thresholds,1);
    concordance = nan(n_thresholds,1);
    ci_lower = nan(n_thresholds,1);
    ci_upper = nan(n_thresholds,1);

    for k = 1:n_thresholds
        threshold = thresholds(k);
        outside = abs(x) > threshold | abs(y) > threshold;
        use = outside & ~baseline;
        x_used = x(use);
        y_used = y(use);
        n = numel(x_used);
        count = sum(x_used.*y_used > 0);

        n_included(k) = n;
        concordant_count(k) = count;
        if n > 0
            p = count/n;
            concordance(k) = p;
            z95 = 1.96;
            denominator = 1 + z95^2/n;
            center = (p + z95^2/(2*n))/denominator;
            half_width = z95*sqrt(p*(1-p)/n + z95^2/(4*n^2))/denominator;
            ci_lower(k) = center-half_width;
            ci_upper(k) = center+half_width;
        end
    end

    results = table( ...
        repmat(string(device),n_thresholds,1), thresholds(:), ...
        concordant_count, n_included, concordance, ci_lower, ci_upper, ...
        'VariableNames', {'Device','ExclusionThreshold','ConcordantCount', ...
        'NIncluded','Concordance','CI95_Lower','CI95_Upper'});
end

function draw_four_quadrant_panel(ax, dx, dy, is_baseline, results, device)
    valid = ~isnan(dx) & ~isnan(dy);
    x = dy(valid);
    y = dx(valid);
    baseline = is_baseline(valid);
    lim = 80;
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');

    h10 = patch(ax, [-10 10 10 -10], [-10 -10 10 10], [1.0 0.4 0.6], ...
        'FaceAlpha', 0.25, 'EdgeColor', [1.0 0.4 0.6], 'LineWidth', 1.6, ...
        'DisplayName', '\pm10% exclusion');
    h05 = patch(ax, [-5 5 5 -5], [-5 -5 5 5], [1.0 0.0 0.2], ...
        'FaceAlpha', 0.35, 'EdgeColor', [1.0 0.0 0.2], 'LineWidth', 1.6, ...
        'DisplayName', '\pm5% exclusion');

    outside10 = abs(x) > 10 | abs(y) > 10;
    scatter(ax, x(~outside10), y(~outside10), 28, [0.4 0.4 0.4], ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.4, ...
        'HandleVisibility', 'off');
    scatter(ax, x(outside10), y(outside10), 28, [0 0 0], ...
        'filled', 'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.6, ...
        'HandleVisibility', 'off');
    scatter(ax, x(baseline), y(baseline), 34, [0.2 0.45 0.8], ...
        'o', 'LineWidth', 1, 'HandleVisibility', 'off');

    plot(ax, [-lim lim], [0 0], 'k-', 'HandleVisibility', 'off');
    plot(ax, [0 0], [-lim lim], 'k-', 'HandleVisibility', 'off');
    plot(ax, [-lim lim], [-lim lim], 'k--', 'HandleVisibility', 'off');
    axis(ax, 'equal');
    axis(ax, [-lim lim -lim lim]);
    xlabel(ax, 'SjvO_2 (relative change %)');
    ylabel(ax, 'rSO_2 (relative change %)');
    title(ax, device);

    s5 = results(results.ExclusionThreshold == 5, :);
    s10 = results(results.ExclusionThreshold == 10, :);
    text(ax, 0.72*lim, -0.30*lim, sprintf('%.1f %%',100*s5.Concordance), ...
        'Color', [1.0 0.0 0.2], 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'right');
    text(ax, 0.72*lim, -0.38*lim, sprintf('%.1f %%',100*s10.Concordance), ...
        'Color', [1.0 0.4 0.6], 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'right');
    legend(ax, [h05 h10], 'Location', 'southeast', 'Box', 'off');
end
