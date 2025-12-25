%% ========================================================================
%  SOBOL VISUALIZATION - FIXED VERSION
%  Resolves fill() dimension mismatch error
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
sum_S1 = results.sum_S;
sum_ST = results.sum_ST;

% Ensure time is column vector
time = time(:);

% Parameter names
param_names = {'Temperature', 'Pressure', 'Flow Rate'};
param_short = {'T', 'P', 'F'};

% Colors for consistency (colorblind-friendly)
colors = [0.90 0.35 0.25;  % Red - Temperature
          0.20 0.45 0.70;  % Blue - Pressure
          0.30 0.70 0.30]; % Green - Flow Rate

fprintf('✓ Data loaded successfully\n');
fprintf('  Time points: %d (from %.0f to %.0f min)\n', length(time), min(time), max(time));
fprintf('  Sample size: N = %d\n\n', results.N);

%% ========================================================================
%% FIGURE 1: Total-Order Indices Evolution (MAIN RESULT)
%% ========================================================================

fprintf('Creating Figure 1: Total-Order Evolution...\n');

figure('Position', [100 100 900 600], 'Color', 'w');

% Plot total-order indices
for i = 1:3
    plot(time, ST(:,i), 'o-', 'LineWidth', 2.5, 'MarkerSize', 6, ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
         'DisplayName', param_names{i});
    hold on;
end

% Formatting
xlabel('Time (min)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Total-Order Sobol Index (S_T)', 'FontSize', 14, 'FontWeight', 'bold');
title('Parameter Importance Evolution (N=32,768)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'east', 'FontSize', 12);
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0 max(time)*1.05]);
ylim([0 1.1]);

% Add phase annotations
y_ann = 1.02;
text(50, y_ann, 'Kinetic-Limited', 'FontSize', 10, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 0.9]);
text(600, y_ann, 'Thermodynamic-Limited', 'FontSize', 10, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'BackgroundColor', [0.9 1 1]);

% Save
saveas(gcf, 'Fig1_TotalOrder_Evolution.png');
saveas(gcf, 'Fig1_TotalOrder_Evolution.fig');
fprintf('✓ Figure 1 saved\n\n');

%% ========================================================================
%% FIGURE 2: First-Order vs Total-Order Comparison (FIXED)
%% ========================================================================

fprintf('Creating Figure 2: First vs Total-Order Comparison...\n');

figure('Position', [150 150 1200 500], 'Color', 'w');

for i = 1:3
    subplot(1,3,i);
    
    % Plot both first-order and total-order
    plot(time, S1(:,i), 's--', 'LineWidth', 2, 'MarkerSize', 7, ...
         'Color', colors(i,:)*0.6, 'MarkerFaceColor', colors(i,:)*0.6, ...
         'DisplayName', 'First-Order (S)');
    hold on;
    plot(time, ST(:,i), 'o-', 'LineWidth', 2.5, 'MarkerSize', 7, ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
         'DisplayName', 'Total-Order (ST)');
    
    % Shaded area between (interaction) - FIXED VERSION
    x_fill = [time; flipud(time)];
    y_fill = [S1(:,i); flipud(ST(:,i))];
    
    % Fill only if no errors
    try
        fill(x_fill, y_fill, colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'DisplayName', 'Interactions', 'HandleVisibility', 'off');
    catch
        % If fill fails, just skip it
        fprintf('  Note: Skipping shaded area for %s (dimension issue)\n', param_names{i});
    end
    
    % Formatting
    xlabel('Time (min)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Sobol Index', 'FontSize', 12, 'FontWeight', 'bold');
    title(param_names{i}, 'FontSize', 14, 'FontWeight', 'bold', 'Color', colors(i,:));
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    box on;
    set(gca, 'FontSize', 11, 'LineWidth', 1.5);
    xlim([0 max(time)*1.05]);
    
    % Add note for Temperature
    if i == 1
        text(max(time)*0.5, max(ST(:,i))*0.7, ...
             {'First-order unreliable', '(parameter redundancy)'}, ...
             'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 0.9]);
    end
end

sgtitle('First-Order vs Total-Order Indices (N=32,768)', 'FontSize', 16, 'FontWeight', 'bold');

% Save
saveas(gcf, 'Fig2_FirstOrder_vs_TotalOrder.png');
saveas(gcf, 'Fig2_FirstOrder_vs_TotalOrder.fig');
fprintf('✓ Figure 2 saved\n\n');

%% ========================================================================
%% FIGURE 3: Parameter Ranking at t=1100
%% ========================================================================

fprintf('Creating Figure 3: Parameter Ranking...\n');

figure('Position', [200 200 800 600], 'Color', 'w');

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
ylabel('Total-Order Index', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Parameter Importance at t = %d min', t_actual), ...
      'FontSize', 16, 'FontWeight', 'bold');
grid on;
box on;
set(gca, 'FontSize', 13, 'LineWidth', 2, 'FontWeight', 'bold');
ylim([0 0.6]);

% Add percentage labels
for i = 1:3
    text(i, ST_sorted(i)+0.025, sprintf('%.1f%%', ST_sorted(i)*100), ...
         'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Save
saveas(gcf, 'Fig3_Parameter_Ranking.png');
saveas(gcf, 'Fig3_Parameter_Ranking.fig');
fprintf('✓ Figure 3 saved\n\n');

%% ========================================================================
%% FIGURE 4: Second and Third-Order Interactions
%% ========================================================================

fprintf('Creating Figure 4: Interactions...\n');

figure('Position', [250 250 1000 600], 'Color', 'w');

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
ylabel('Second-Order Indices', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0 0.10]);
set(gca, 'YColor', 'k');

% Plot third-order
yyaxis right
plot(time, abs(S3), 'd-', 'LineWidth', 3, 'MarkerSize', 8, ...
     'Color', [0.7 0.0 0.0], 'MarkerFaceColor', [0.7 0.0 0.0], ...
     'DisplayName', '|S3| (T × P × F)');
ylabel('Third-Order Index (Magnitude)', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0 0.30]);
set(gca, 'YColor', [0.7 0.0 0.0]);

% Formatting
xlabel('Time (min)', 'FontSize', 14, 'FontWeight', 'bold');
title('Parameter Interaction Effects', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 12);
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 2);
xlim([0 max(time)*1.05]);

% Save
saveas(gcf, 'Fig4_Interactions.png');
saveas(gcf, 'Fig4_Interactions.fig');
fprintf('✓ Figure 4 saved\n\n');

%% ========================================================================
%% FIGURE 5: Sum of Indices (Quality Check)
%% ========================================================================

fprintf('Creating Figure 5: Quality Check...\n');

figure('Position', [300 300 900 600], 'Color', 'w');

% Plot sums
plot(time, sum_S1, 's-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
     'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10], ...
     'DisplayName', 'Sum(S_i) - First-Order');
hold on;
plot(time, sum_ST, 'o-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
     'Color', [0.00 0.45 0.74], 'MarkerFaceColor', [0.00 0.45 0.74], ...
     'DisplayName', 'Sum(ST_i) - Total-Order');

% Reference line at 1.0
yline(1.0, 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (1.0)');

% Shaded acceptable region
fill([0 max(time) max(time) 0], [0.95 0.95 1.05 1.05], ...
     [0.9 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
     'DisplayName', '±5% Tolerance', 'HandleVisibility', 'off');

% Formatting
xlabel('Time (min)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Sum of Indices', 'FontSize', 14, 'FontWeight', 'bold');
title('Variance Decomposition Quality (N=32,768)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 12);
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0 max(time)*1.05]);
ylim([0.9 1.3]);

% Save
saveas(gcf, 'Fig5_Sum_of_Indices.png');
saveas(gcf, 'Fig5_Sum_of_Indices.fig');
fprintf('✓ Figure 5 saved\n\n');

%% ========================================================================
%% Summary
%% ========================================================================

fprintf('========================================================================\n');
fprintf('VISUALIZATION COMPLETE!\n');
fprintf('========================================================================\n');
fprintf('\nFigures created:\n');
fprintf('  1. Fig1_TotalOrder_Evolution.png/.fig\n');
fprintf('  2. Fig2_FirstOrder_vs_TotalOrder.png/.fig\n');
fprintf('  3. Fig3_Parameter_Ranking.png/.fig\n');
fprintf('  4. Fig4_Interactions.png/.fig\n');
fprintf('  5. Fig5_Sum_of_Indices.png/.fig\n');
fprintf('\n');

% Print summary table
fprintf('KEY RESULTS (t = 1100 min):\n');
fprintf('----------------------------\n');
idx_1100 = find(time == 1100, 1);
if isempty(idx_1100)
    [~, idx_1100] = min(abs(time - 1100));
end

fprintf('Total-Order Indices:\n');
fprintf('  Temperature: %.4f (%.1f%%)\n', ST(idx_1100,1), ST(idx_1100,1)*100);
fprintf('  Pressure:    %.4f (%.1f%%)\n', ST(idx_1100,2), ST(idx_1100,2)*100);
fprintf('  Flow Rate:   %.4f (%.1f%%)\n', ST(idx_1100,3), ST(idx_1100,3)*100);
fprintf('  Sum:         %.4f\n', sum_ST(idx_1100));
fprintf('\n');

fprintf('Parameter Ranking:\n');
[~, rank] = sort(ST(idx_1100,:), 'descend');
for i = 1:3
    fprintf('  #%d: %s\n', i, param_names{rank(i)});
end
fprintf('\n');

fprintf('========================================================================\n');
fprintf('Ready for thesis!\n');
fprintf('========================================================================\n');