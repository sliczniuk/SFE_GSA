%% ========================================================================
%  SOBOL SENSITIVITY ANALYSIS - ESSENTIAL FIGURES ONLY
%  N = 32,768 Results - Quick Visualization
%  
%  Creates 4 key figures for thesis
%  ========================================================================

clear; close all; clc;

%% Load Data
fprintf('Loading N=32,768 results...\n');
load('GSA_results_N32768_2025-12-23_06-18.mat');

% Extract data
time = results.time_points;
S1 = results.first_order;
ST = results.total_order;
S2 = results.second_order;
S3 = results.third_order;
sum_ST = results.sum_ST;

% Parameter names
param_names = {'Temperature', 'Pressure', 'Flow Rate'};

% Colors (colorblind-friendly)
colors = [0.90 0.35 0.25;  % Red - Temperature
          0.20 0.45 0.70;  % Blue - Pressure
          0.30 0.70 0.30]; % Green - Flow Rate

%% ========================================================================
%% FIGURE 1: Total-Order Evolution (MAIN THESIS FIGURE)
%% ========================================================================

figure('Position', [100 100 900 650], 'Color', 'w');

% Plot total-order indices
for i = 1:3
    plot(time, ST(:,i), 'o-', 'LineWidth', 3, 'MarkerSize', 8, ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
         'DisplayName', param_names{i});
    hold on;
end

% Formatting
xlabel('Extraction Time (min)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Total-Order Sensitivity Index (S_T_i)', 'FontSize', 16, 'FontWeight', 'bold');
title('Parameter Importance Evolution', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'east', 'FontSize', 14, 'Box', 'on');
grid on;
box on;
set(gca, 'FontSize', 14, 'LineWidth', 2);
xlim([0 max(time)+50]);
ylim([0 1.05]);

% Add phase annotations
annotation('textbox', [0.15 0.88 0.2 0.08], 'String', 'Kinetic Phase', ...
    'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', [0.8 0.8 0], ...
    'BackgroundColor', [1 1 0.9], 'HorizontalAlignment', 'center');
annotation('textbox', [0.65 0.88 0.25 0.08], 'String', 'Thermodynamic Phase', ...
    'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', [0 0.6 0.8], ...
    'BackgroundColor', [0.9 1 1], 'HorizontalAlignment', 'center');

% Save high-resolution
print('Fig_TotalOrder_Main', '-dpng', '-r300');
saveas(gcf, 'Fig_TotalOrder_Main.fig');
fprintf('✓ Figure 1 saved: Total-Order Evolution (MAIN)\n');

%% ========================================================================
%% FIGURE 2: Parameter Ranking at t=1100 (SUMMARY FIGURE)
%% ========================================================================

figure('Position', [150 150 800 650], 'Color', 'w');

% Find t=1100
[~, idx] = min(abs(time - 1100));
t_actual = time(idx);

% Data
ST_values = ST(idx,:);
[ST_sorted, sort_idx] = sort(ST_values, 'descend');
params_sorted = param_names(sort_idx);
colors_sorted = colors(sort_idx,:);

% Create bar chart
b = bar(1:3, ST_sorted, 'FaceColor', 'flat');
b.CData = colors_sorted;
b.LineWidth = 1.5;

% Formatting
xticks(1:3);
xticklabels(params_sorted);
ylabel('Total-Order Index', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Parameter Importance at t = %d min', t_actual), ...
      'FontSize', 18, 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'FontWeight', 'bold');
ylim([0 0.6]);

% Add percentage labels
for i = 1:3
    text(i, ST_sorted(i)+0.025, sprintf('%.1f%%', ST_sorted(i)*100), ...
         'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Add ranking annotation
annotation('textbox', [0.65 0.75 0.25 0.15], ...
    'String', {sprintf('Ranking:'), sprintf('1. %s', params_sorted{1}), ...
               sprintf('2. %s', params_sorted{2}), sprintf('3. %s', params_sorted{3})}, ...
    'FontSize', 13, 'FontWeight', 'bold', 'EdgeColor', [0 0.5 0], ...
    'BackgroundColor', [0.9 1 0.9], 'LineWidth', 2);

% Save
print('Fig_Ranking_Summary', '-dpng', '-r300');
saveas(gcf, 'Fig_Ranking_Summary.fig');
fprintf('✓ Figure 2 saved: Parameter Ranking at t=%d\n', t_actual);

%% ========================================================================
%% FIGURE 3: Interaction Effects (SECOND & THIRD ORDER)
%% ========================================================================

figure('Position', [200 200 1000 650], 'Color', 'w');

% Extract interactions
S2_TP = squeeze(S2(:,1,2));
S2_PF = squeeze(S2(:,2,3));

% Plot second-order
yyaxis left
plot(time, S2_TP, 'o-', 'LineWidth', 2.5, 'MarkerSize', 7, ...
     'Color', [0.5 0.0 0.5], 'MarkerFaceColor', [0.5 0.0 0.5], ...
     'DisplayName', 'S2_{TP} (T × P)');
hold on;
plot(time, S2_PF, 's-', 'LineWidth', 2.5, 'MarkerSize', 7, ...
     'Color', [0.0 0.6 0.6], 'MarkerFaceColor', [0.0 0.6 0.6], ...
     'DisplayName', 'S2_{PF} (P × F)');
ylabel('Second-Order Indices', 'FontSize', 15, 'FontWeight', 'bold');
ylim([0 0.10]);
set(gca, 'YColor', 'k');

% Plot third-order
yyaxis right
plot(time, abs(S3), 'd-', 'LineWidth', 3, 'MarkerSize', 8, ...
     'Color', [0.7 0.0 0.0], 'MarkerFaceColor', [0.7 0.0 0.0], ...
     'DisplayName', '|S3| (T × P × F)');
ylabel('Third-Order Index (Magnitude)', 'FontSize', 15, 'FontWeight', 'bold');
ylim([0 0.30]);
set(gca, 'YColor', [0.7 0.0 0.0]);

% Formatting
xlabel('Extraction Time (min)', 'FontSize', 16, 'FontWeight', 'bold');
title('Parameter Interaction Effects', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 13);
grid on;
box on;
set(gca, 'FontSize', 14, 'LineWidth', 2);
xlim([0 max(time)+50]);

% Add annotation for negative S3
annotation('textbox', [0.55 0.25 0.35 0.15], ...
    'String', {'S3 is NEGATIVE:', 'Three-way compensation', ...
               '(Joint < Sum of effects)'}, ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.7 0.0 0.0], ...
    'EdgeColor', [0.7 0.0 0.0], 'BackgroundColor', [1 0.95 0.95], 'LineWidth', 2);

% Save
print('Fig_Interactions', '-dpng', '-r300');
saveas(gcf, 'Fig_Interactions.fig');
fprintf('✓ Figure 3 saved: Interaction Effects\n');

%% ========================================================================
%% FIGURE 4: Variance Budget at Key Times (4-panel)
%% ========================================================================

figure('Position', [250 250 1200 800], 'Color', 'w');

time_points = [300, 600, 1100, 1500];

for i = 1:4
    [~, idx] = min(abs(time - time_points(i)));
    
    subplot(2, 2, i);
    
    % Main effects
    main_effects = ST(idx,:);
    
    % Interaction contributions
    interact_T = max(0, ST(idx,1) - S1(idx,1));
    interact_P = max(0, ST(idx,2) - S1(idx,2));
    interact_F = max(0, ST(idx,3) - S1(idx,3));
    interactions = [interact_T, interact_P, interact_F];
    
    % Stack plot
    h = bar(1:3, [S1(idx,:); interactions]', 'stacked');
    h(1).FaceColor = 'flat';
    h(1).CData = colors;
    h(2).FaceColor = [0.7 0.7 0.7];
    h(2).FaceAlpha = 0.5;
    
    % Formatting
    xticks(1:3);
    xticklabels({'T', 'P', 'F'});
    ylabel('Sensitivity Index', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('t = %d min', time(idx)), 'FontSize', 15, 'FontWeight', 'bold');
    legend('Main Effect', 'Interactions', 'Location', 'northwest', 'FontSize', 10);
    grid on;
    box on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);
    ylim([0 0.8]);
    
    % Add total values
    for j = 1:3
        text(j, ST(idx,j)+0.03, sprintf('%.2f', ST(idx,j)), ...
             'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
end

sgtitle('Variance Decomposition Over Time', 'FontSize', 18, 'FontWeight', 'bold');

% Save
print('Fig_Variance_Budget', '-dpng', '-r300');
saveas(gcf, 'Fig_Variance_Budget.fig');
fprintf('✓ Figure 4 saved: Variance Budget\n');

%% ========================================================================
%% Print Summary Table
%% ========================================================================

fprintf('\n');
fprintf('========================================================================\n');
fprintf('SUMMARY TABLE FOR THESIS (t = 1100 min)\n');
fprintf('========================================================================\n');
fprintf('Parameter      | Total-Order | Percentage | Ranking\n');
fprintf('------------------------------------------------------------------------\n');

[~, idx_1100] = min(abs(time - 1100));
ST_1100 = ST(idx_1100,:);
[ST_sorted, rank] = sort(ST_1100, 'descend');

for i = 1:3
    idx_param = rank(i);
    fprintf('%-14s | %8.4f    | %6.1f%%    | #%d\n', ...
            param_names{idx_param}, ST_1100(idx_param), ...
            ST_1100(idx_param)*100, i);
end

fprintf('------------------------------------------------------------------------\n');
fprintf('Sum(ST)        | %8.4f    | %6.1f%%    |\n', ...
        sum_ST(idx_1100), sum_ST(idx_1100)*100);
fprintf('========================================================================\n');

fprintf('\nSecond-Order Interactions:\n');
fprintf('  T × P:  %.4f (%.1f%%)\n', S2(idx_1100,1,2), S2(idx_1100,1,2)*100);
fprintf('  P × F:  %.4f (%.1f%%)\n', S2(idx_1100,2,3), S2(idx_1100,2,3)*100);
fprintf('  T × F:  %.4f (%.1f%%)\n', S2(idx_1100,1,3), S2(idx_1100,1,3)*100);
fprintf('\n');

fprintf('Third-Order Interaction:\n');
fprintf('  T × P × F:  %.4f (%.1f%% compensation)\n', S3(idx_1100), abs(S3(idx_1100))*100);
fprintf('\n');

fprintf('========================================================================\n');
fprintf('Essential figures created successfully!\n');
fprintf('Files saved:\n');
fprintf('  - Fig_TotalOrder_Main.png/.fig (Use in RESULTS)\n');
fprintf('  - Fig_Ranking_Summary.png/.fig (Use in ABSTRACT/SUMMARY)\n');
fprintf('  - Fig_Interactions.png/.fig (Use in DISCUSSION)\n');
fprintf('  - Fig_Variance_Budget.png/.fig (Use in ANALYSIS)\n');
fprintf('========================================================================\n');

%% Generate LaTeX table code
fid = fopen('Table_Summary.tex', 'w');
fprintf(fid, '%% LaTeX table for thesis\n');
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Parameter Sensitivity at t = %d min (N = 32,768)}\n', time(idx_1100));
fprintf(fid, '\\label{tab:sensitivity}\n');
fprintf(fid, '\\begin{tabular}{lccc}\n');
fprintf(fid, '\\hline\\hline\n');
fprintf(fid, 'Parameter & Total-Order Index & Percentage & Rank \\\\\n');
fprintf(fid, '\\hline\n');

for i = 1:3
    idx_param = rank(i);
    fprintf(fid, '%s & %.4f & %.1f\\%% & %d \\\\\n', ...
            param_names{idx_param}, ST_1100(idx_param), ...
            ST_1100(idx_param)*100, i);
end

fprintf(fid, '\\hline\n');
fprintf(fid, 'Sum & %.4f & %.1f\\%% & --- \\\\\n', ...
        sum_ST(idx_1100), sum_ST(idx_1100)*100);
fprintf(fid, '\\hline\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);

fprintf('\n✓ LaTeX table saved: Table_Summary.tex\n');

fprintf('\n========================================================================\n');
fprintf('QUICK START COMPLETE!\n');
fprintf('Ready for thesis preparation.\n');
fprintf('========================================================================\n');