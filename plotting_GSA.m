% Visualize Sobol Sensitivity Analysis Results
% This script creates comprehensive visualizations of sensitivity indices evolution

clear; close all; clc;

%% Data Entry
load GSA_results_N32768_2025-12-23_06-18

time = results.time_points';

S1_T = results.first_order(:,1);
S1_P = results.first_order(:,2);
S1_F = results.first_order(:,3);

ST_T = results.total_order(:,1);
ST_P = results.total_order(:,2);
ST_F = results.total_order(:,3);

mean_yield = results.output_mean;
std_yield  = results.output_std;

sum_S1     = results.sum_S;
sum_ST     = results.sum_ST;

%%
%load GSA_results_N10000_2025-12-12_16-06

%%
%{
% Time points [min]
time = [5, 15, 30, 60, 90, 120, 150, 240, 300, 450, 600, 750, 900, 1200, 1500, 2000]';

% First-order Sobol indices
S1_T = [-0.0214, -0.0276, -0.0282, -0.0095, 0.0076, 0.0203, 0.0303, ...
        0.0521, 0.0635, 0.0884, 0.1123, 0.1369, 0.1630, 0.2208, ...
        0.2893, 0.4398]';  % Temperature
S1_P = [0.0742, 0.0224, -0.0276, 0.0085, 0.0678, 0.1162, 0.1546, ...
        0.2354, 0.2743, 0.3514, 0.4174, 0.4808, 0.5441, 0.6712, ...
        0.7937, 0.9533]';  % Pressure
S1_F = [0.8467, 0.9413, 0.9995, 0.9615, 0.9003, 0.8534, 0.8189, ...
        0.7579, 0.7360, 0.7092, 0.7000, 0.6972, 0.6961, 0.6895, ...
        0.6687, 0.5832]';  % Flow rate

%S1_T = max(0,S1_T);
%S1_F = max(0,S1_F);
%S1_P = max(0,S1_P);

% Total-order Sobol indices
ST_T = [0.0195, 0.0068, 0.0016, 0.0067, 0.0137, 0.0191, 0.0232, ...
        0.0314, 0.0352, 0.0425, 0.0490, 0.0557, 0.0627, 0.0780, ...
        0.0949, 0.1256]';  % Temperature
ST_P = [0.1355, 0.0553, 0.0137, 0.0561, 0.1152, 0.1608, 0.1953, ...
        0.2614, 0.2893, 0.3361, 0.3690, 0.3965, 0.4219, 0.4700, ...
        0.5165, 0.5899]';  % Pressure
ST_F = [0.9064, 0.9678, 1.0282, 0.9883, 0.9227, 0.8702, 0.8294, ...
        0.7480, 0.7122, 0.6505, 0.6082, 0.5754, 0.5482, 0.5045, ...
        0.4706, 0.4290]';  % Flow rate

% Output statistics
mean_yield = [0.0213, 0.3364, 0.8522, 1.3949, 1.6970, 1.9004, 2.0502, ...
              2.3378, 2.4585, 2.6475, 2.7553, 2.8234, 2.8689, 2.9235, ...
              2.9529, 2.9773]';
std_yield = [0.0125, 0.1262, 0.1667, 0.1794, 0.1874, 0.1911, 0.1918, ...
             0.1849, 0.1767, 0.1530, 0.1299, 0.1096, 0.0924, 0.0660, ...
             0.0477, 0.0286]';

% Sum of indices
sum_S1 = [0.8995, 0.9361, 0.9437, 0.9605, 0.9757, 0.9900, 1.0039, ...
          1.0454, 1.0738, 1.1489, 1.2297, 1.3150, 1.4031, 1.5815, ...
          1.7517, 1.9763]';
sum_ST = [1.0613, 1.0299, 1.0435, 1.0512, 1.0516, 1.0500, 1.0478, ...
          1.0408, 1.0366, 1.0292, 1.0262, 1.0276, 1.0327, 1.0526, ...
          1.0820, 1.1444]';
%}

%% Calculate confidence bounds
% Using ±1σ (68% confidence)
upper_1sigma = mean_yield + std_yield;
lower_1sigma = mean_yield - std_yield;

% Using ±2σ (95% confidence)
upper_2sigma = mean_yield + 2*std_yield;
lower_2sigma = mean_yield - 2*std_yield;

% Using ±3σ (99.7% confidence)
upper_3sigma = mean_yield + 3*std_yield;
lower_3sigma = mean_yield - 3*std_yield;

%% Figure 1: First-order Sobol indices evolution
figure('Position', [100, 100, 1200, 500]);

subplot(1,2,1)
plot(time, S1_T, 'r', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Temperature'); hold on;
plot(time, S1_P, 'b', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Pressure');
plot(time, S1_F, 'g', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Flow rate');
plot(time, zeros(size(time)), 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel('Time [min]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('First-order Sobol index, $S_i$', 'FontSize', 14, 'FontWeight', 'bold');
%title('First-order Sensitivity Indices', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 14); legend box off
grid on;
set(gca, 'FontSize', 16);
xlim([0, max(time)]);

subplot(1,2,2)

plot(time, ST_T, 'r', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Temperature'); hold on;
plot(time, ST_P, 'b', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Pressure');
plot(time, ST_F, 'g', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Flow rate');
xlabel('Time [min]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Total-order Sobol index, $S_{Ti}$', 'FontSize', 14, 'FontWeight', 'bold');
%title('Evolution of Total-order Sensitivity Indices', 'FontSize', 10, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 12); legend box off;
grid on;
set(gca, 'FontSize', 16);
xlim([0, max(time)]);
exportgraphics(figure(1),['Sobol.png'], "Resolution",500); close all

%% Figure 2: Multiple confidence levels
figure('Position', [120, 120, 1000, 600]);

% Plot 3σ (outermost, lightest)
fill([time; flipud(time)], [upper_3sigma; flipud(lower_3sigma)], ...
    [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;

% Plot 2σ (middle)
fill([time; flipud(time)], [upper_2sigma; flipud(lower_2sigma)], ...
    [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot 1σ (innermost, darkest)
fill([time; flipud(time)], [upper_1sigma; flipud(lower_1sigma)], ...
    [0.3 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

% Plot mean yield curve
yyaxis left
plot(time, mean_yield, 'k-', 'LineWidth', 3, 'Color', [0 0 0]);
plot(time, mean_yield, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [0 0 0], ...
    'MarkerEdgeColor', [0 0 0])

% Formatting
xlabel('Time [min]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Yield [g]', 'FontSize', 14, 'FontWeight', 'bold');

xlim([0, max(time)*1.02]);
ylim([0, max(upper_3sigma)*1.1]);

yyaxis right
plot(time, std_yield, 'b', 'LineWidth', 3)
ylabel('$\sigma$', 'FontSize', 14, 'FontWeight', 'bold');

%title('Extraction Yield with Multiple Confidence Levels', 'FontSize', 16, 'FontWeight', 'bold');
legend('99.7\% CI', '95\% CI', '68\% CI', ...
    'Mean yield [g]', '','$\sigma$' ,'Location', 'east', 'FontSize', 14); legend box off;
grid on;
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
exportgraphics(figure(1),['CI.png'], "Resolution",500); close all

%% Figure 4: Interaction effects
figure('Position', [160, 160, 1000, 600]);

interaction = ST_T - S1_T;  % Temperature interactions
interaction_P = ST_P - S1_P;  % Pressure interactions
interaction_F = ST_F - S1_F;  % Flow rate interactions

plot(time, interaction, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Temperature'); hold on;
plot(time, interaction_P, 'b-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Pressure');
plot(time, interaction_F, 'g-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Flow rate');
plot(time, zeros(size(time)), 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel('Time [min]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Interaction index ($S_{Ti} - S_i$)', 'FontSize', 14, 'FontWeight', 'bold');
title('Parameter Interaction Effects', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14); legend Box off;
grid on;
set(gca, 'FontSize', 16);
xlim([0, max(time)]);
exportgraphics(figure(1),['Interaction_effects.png'], "Resolution",500); close all

%% Save figures
% Uncomment to save figures
% print(1, 'sobol_first_order', '-dpng', '-r300');
% print(2, 'sobol_total_order', '-dpng', '-r300');
% print(3, 'sobol_stacked_area', '-dpng', '-r300');
% print(4, 'sobol_interactions', '-dpng', '-r300');
% print(5, 'sobol_sum_indices', '-dpng', '-r300');
% print(6, 'sobol_output_stats', '-dpng', '-r300');
% print(7, 'sobol_key_timepoints', '-dpng', '-r300');
% print(8, 'sobol_dashboard', '-dpng', '-r300');

fprintf('Visualization complete!\n');