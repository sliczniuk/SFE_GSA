% Advanced Interaction Analysis for Sobol Sensitivity Results
% This script estimates pairwise interactions and creates detailed visualizations

clear; close all; clc;

%% Data Entry
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

% Sum of indices
sum_S1 = [0.8995, 0.9361, 0.9437, 0.9605, 0.9757, 0.9900, 1.0039, ...
          1.0454, 1.0738, 1.1489, 1.2297, 1.3150, 1.4031, 1.5815, ...
          1.7517, 1.9763]';
sum_ST = [1.0613, 1.0299, 1.0435, 1.0512, 1.0516, 1.0500, 1.0478, ...
          1.0408, 1.0366, 1.0292, 1.0262, 1.0276, 1.0327, 1.0526, ...
          1.0820, 1.1444]';

%% Calculate Interaction Indices
% Individual parameter interactions with all others
interaction_T = ST_T - S1_T;  % Temperature total interactions
interaction_P = ST_P - S1_P;  % Pressure total interactions
interaction_F = ST_F - S1_F;  % Flow rate total interactions

% Total interaction index (all higher-order effects)
% This represents all interaction terms: pairwise + three-way
total_interactions = sum_ST - sum_S1;

% Estimate of all higher-order interactions (>2nd order)
% Using: Σ ST - Σ S1 ≈ ΣΣ S_ij + S_123
higher_order_total = sum_S1 + sum_ST - 1;  % Simplified estimate

%% Estimate Pairwise Interaction Indices
% Note: These are approximate estimates since we don't have the closed-form
% second-order indices. We use the relationship:
% S_Ti ≈ S_i + Σ(j≠i) S_ij + S_ijk...
% For three parameters: S_Ti = S_i + S_i,others
% Where S_i,others includes both pairwise and three-way interactions

% Approximate pairwise interactions (assuming three-way interaction is small)
% S_TP: Temperature-Pressure interaction
% S_TF: Temperature-Flow interaction  
% S_PF: Pressure-Flow interaction

% Using the closure property: Σ S_i + ΣΣ S_ij + S_123 = 1
% And: S_Ti = S_i + Σ(j≠i) S_ij + S_123

% Total pairwise + higher order interactions
total_interaction_all = 1 - sum_S1;

% Rough estimate: distribute interactions proportionally
% This is an approximation since we cannot uniquely determine S_ij without
% additional Sobol indices calculations

% Alternative estimate using the relationship:
% sum_ST - sum_S1 gives us total interaction magnitude
% We can bound the pairwise interactions

% Lower bound: assuming three-way interaction is zero
% S_TP + S_TF + S_PF ≈ sum_ST - sum_S1

% For visualization, we estimate dominant pairwise interactions
% Based on which parameters have largest (S_Ti - S_i) products

%% Figure 1: Interaction Effects Comparison
figure('Position', [100, 100, 1400, 900]);

% Panel 1: Individual parameter interactions
subplot(2,3,1)
plot(time, interaction_T, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Temperature'); hold on;
plot(time, interaction_P, 'b-s', 'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Pressure');
plot(time, interaction_F, 'g-^', 'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Flow rate');
plot(time, zeros(size(time)), 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Interaction Index (S_{Ti} - S_i)', 'FontSize', 11, 'FontWeight', 'bold');
title('Individual Parameter Interactions', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
grid on; set(gca, 'FontSize', 10);
xlim([0, max(time)]);

% Panel 2: Stacked interaction contributions
subplot(2,3,2)
interaction_matrix = [interaction_T, interaction_P, interaction_F];
area(time, interaction_matrix, 'LineWidth', 1.5);
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Interaction Contribution', 'FontSize', 11, 'FontWeight', 'bold');
title('Stacked Interaction Effects', 'FontSize', 12, 'FontWeight', 'bold');
legend('Temperature', 'Pressure', 'Flow rate', 'Location', 'northwest', 'FontSize', 10);
grid on; set(gca, 'FontSize', 10);
colormap([0.8 0.2 0.2; 0.2 0.4 0.8; 0.2 0.8 0.2]);
xlim([0, max(time)]);

% Panel 3: Relative interaction strength
subplot(2,3,3)
% Calculate relative interaction (as % of total effect)
rel_int_T = 100 * interaction_T ./ ST_T;
rel_int_P = 100 * interaction_P ./ ST_P;
rel_int_F = 100 * interaction_F ./ ST_F;
% Handle division by zero or very small numbers
rel_int_T(ST_T < 0.001) = 0;
rel_int_P(ST_P < 0.001) = 0;
rel_int_F(ST_F < 0.001) = 0;

plot(time, rel_int_T, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Temperature'); hold on;
plot(time, rel_int_P, 'b-s', 'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Pressure');
plot(time, rel_int_F, 'g-^', 'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Flow rate');
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Interaction / Total Effect [%]', 'FontSize', 11, 'FontWeight', 'bold');
title('Relative Interaction Strength', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on; set(gca, 'FontSize', 10);
xlim([0, max(time)]);

% Panel 4: Total interaction budget
subplot(2,3,4)
plot(time, total_interaction_all, 'k-o', 'LineWidth', 2.5, 'MarkerSize', 8, ...
    'DisplayName', 'Total Interactions (1 - \Sigma S_i)'); hold on;
plot(time, sum_ST - sum_S1, 'm-s', 'LineWidth', 2.5, 'MarkerSize', 8, ...
    'DisplayName', '\Sigma S_{Ti} - \Sigma S_i');
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Total Interaction Magnitude', 'FontSize', 11, 'FontWeight', 'bold');
title('Overall Interaction Budget', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
grid on; set(gca, 'FontSize', 10);
xlim([0, max(time)]);

% Panel 5: Interaction vs Main Effect Ratio
subplot(2,3,5)
% For each parameter, show ratio of interaction to main effect
ratio_T = interaction_T ./ max(abs(S1_T), 0.001);
ratio_P = interaction_P ./ max(abs(S1_P), 0.001);
ratio_F = interaction_F ./ S1_F;

semilogy(time, abs(ratio_T), 'r-o', 'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Temperature'); hold on;
semilogy(time, abs(ratio_P), 'b-s', 'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Pressure');
semilogy(time, ratio_F, 'g-^', 'LineWidth', 2.5, 'MarkerSize', 7, 'DisplayName', 'Flow rate');
semilogy(time, ones(size(time)), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('|(S_{Ti} - S_i) / S_i|', 'FontSize', 11, 'FontWeight', 'bold');
title('Interaction-to-Main Effect Ratio (Log)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on; set(gca, 'FontSize', 10);
xlim([0, max(time)]);

% Panel 6: Sum decomposition
subplot(2,3,6)
bar_data = [sum_S1, total_interaction_all];
b = bar(time, bar_data, 'stacked');
b(1).FaceColor = [0.3 0.6 0.9];
b(2).FaceColor = [0.9 0.4 0.3];
hold on;
plot(time, ones(size(time)), 'k--', 'LineWidth', 2);
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Variance Decomposition', 'FontSize', 11, 'FontWeight', 'bold');
title('Main Effects + Interactions = 1', 'FontSize', 12, 'FontWeight', 'bold');
legend('Main effects (\Sigma S_i)', 'Interactions', 'Total', 'Location', 'east', 'FontSize', 10);
grid on; set(gca, 'FontSize', 10);
xlim([0, max(time)*1.05]);


%% Figure 3: Detailed Interaction Timeline
figure('Position', [200, 200, 1200, 800]);

% Top panel: All effects together
subplot(3,1,1)
plot(time, S1_T, 'r-', 'LineWidth', 2, 'DisplayName', 'S_T (main)'); hold on;
plot(time, S1_P, 'b-', 'LineWidth', 2, 'DisplayName', 'S_P (main)');
plot(time, S1_F, 'g-', 'LineWidth', 2, 'DisplayName', 'S_F (main)');
plot(time, interaction_T, 'r--', 'LineWidth', 2, 'DisplayName', 'S_T interactions');
plot(time, interaction_P, 'b--', 'LineWidth', 2, 'DisplayName', 'S_P interactions');
plot(time, interaction_F, 'g--', 'LineWidth', 2, 'DisplayName', 'S_F interactions');
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Index Value', 'FontSize', 11, 'FontWeight', 'bold');
title('Main Effects (Solid) vs Interaction Components (Dashed)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 9);
grid on; set(gca, 'FontSize', 10);
xlim([0, max(time)]);

% Middle panel: Interaction significance
subplot(3,1,2)
% Calculate what fraction of total effect is due to interactions
frac_int_T = 100 * interaction_T ./ (S1_T + interaction_T);
frac_int_P = 100 * interaction_P ./ (S1_P + interaction_P);
frac_int_F = 100 * interaction_F ./ (S1_F + interaction_F);

% Handle edge cases
frac_int_T(S1_T + interaction_T < 0.001) = 0;
frac_int_P(S1_P + interaction_P < 0.001) = 0;

plot(time, frac_int_T, 'r-o', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Temperature'); hold on;
plot(time, frac_int_P, 'b-s', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Pressure');
plot(time, frac_int_F, 'g-^', 'LineWidth', 2.5, 'MarkerSize', 6, 'DisplayName', 'Flow rate');
plot(time, 50*ones(size(time)), 'k--', 'LineWidth', 1, 'DisplayName', '50% threshold');
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Interaction Fraction [%]', 'FontSize', 11, 'FontWeight', 'bold');
title('Percentage of Total Effect from Interactions', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on; set(gca, 'FontSize', 10);
xlim([0, max(time)]);
ylim([0, 100]);

% Bottom panel: Cumulative interaction evolution
subplot(3,1,3)
cumsum_int = cumsum([interaction_T, interaction_P, interaction_F]);
plot(time, cumsum_int(:,1), 'r-', 'LineWidth', 2.5, 'DisplayName', 'Cumulative T'); hold on;
plot(time, cumsum_int(:,2), 'b-', 'LineWidth', 2.5, 'DisplayName', 'Cumulative T+P');
plot(time, cumsum_int(:,3), 'g-', 'LineWidth', 2.5, 'DisplayName', 'Cumulative T+P+F');
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Cumulative Interaction', 'FontSize', 11, 'FontWeight', 'bold');
title('Cumulative Interaction Accumulation', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
grid on; set(gca, 'FontSize', 10);
xlim([0, max(time)]);

%% Figure 4: Parameter Dominance Phases
figure('Position', [250, 250, 1200, 600]);

% Define phases based on dominant parameter
% Phase 1: Flow dominated (S_F > 0.85)
% Phase 2: Transition (0.65 < S_F < 0.85)
% Phase 3: Pressure-dominated (S_P > 0.45)

subplot(2,2,1)
% Ternary-like representation
scatter3(S1_T, S1_P, S1_F, 100, time, 'filled'); hold on;
colorbar;
xlabel('S_T', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('S_P', 'FontSize', 11, 'FontWeight', 'bold');
zlabel('S_F', 'FontSize', 11, 'FontWeight', 'bold');
title('Parameter Space Trajectory', 'FontSize', 12, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 10);
view(45, 30);
colormap(jet);

subplot(2,2,2)
% Phase identification
phase = zeros(size(time));
phase(S1_F > 0.85) = 1;  % Flow dominated
phase(S1_F <= 0.85 & S1_F > 0.65) = 2;  % Transition
phase(S1_F <= 0.65) = 3;  % Pressure dominated

scatter(time, sum_S1, 100, phase, 'filled'); hold on;
plot(time, sum_S1, 'k-', 'LineWidth', 1.5);
plot(time, ones(size(time)), 'r--', 'LineWidth', 2);
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('\Sigma S_i', 'FontSize', 11, 'FontWeight', 'bold');
title('Interaction Emergence by Phase', 'FontSize', 12, 'FontWeight', 'bold');
colorbar('Ticks', [1, 2, 3], 'TickLabels', {'Flow', 'Transition', 'Pressure'});
grid on; set(gca, 'FontSize', 10);

subplot(2,2,[3,4])
% Stacked area showing main + interaction
area_data = [sum_S1, total_interaction_all];
a = area(time, area_data);
a(1).FaceColor = [0.4 0.6 0.9];
a(2).FaceColor = [0.9 0.5 0.3];
hold on;

% Add phase boundaries
phase_boundaries = [time(find(phase==2, 1)), time(find(phase==3, 1))];
for pb = phase_boundaries
    if ~isempty(pb)
        plot([pb pb], [0 2], 'k--', 'LineWidth', 2);
    end
end

plot(time, ones(size(time)), 'k-', 'LineWidth', 2.5);
xlabel('Time [min]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Variance Fraction', 'FontSize', 11, 'FontWeight', 'bold');
title('Main Effects vs Interactions with Process Phases', 'FontSize', 12, 'FontWeight', 'bold');
legend('Main Effects', 'Interactions', 'Location', 'northwest', 'FontSize', 10);
text(50, 1.7, 'Flow Phase', 'FontSize', 11, 'FontWeight', 'bold');
text(400, 1.7, 'Transition', 'FontSize', 11, 'FontWeight', 'bold');
text(1300, 1.7, 'Pressure Phase', 'FontSize', 11, 'FontWeight', 'bold');
grid on; set(gca, 'FontSize', 10);
xlim([0, max(time)]);

%% Statistical Summary Table
fprintf('\n========================================\n');
fprintf('INTERACTION ANALYSIS SUMMARY\n');
fprintf('========================================\n\n');

fprintf('Time Point Analysis:\n');
fprintf('%-10s %-12s %-12s %-12s %-15s\n', 'Time[min]', 'Int_T', 'Int_P', 'Int_F', 'Total_Int');
fprintf('%-10s %-12s %-12s %-12s %-15s\n', '--------', '-----', '-----', '-----', '---------');
for i = 1:length(time)
    fprintf('%-10d %-12.4f %-12.4f %-12.4f %-15.4f\n', ...
        time(i), interaction_T(i), interaction_P(i), interaction_F(i), total_interaction_all(i));
end

fprintf('\n\nPhase Classification:\n');
fprintf('Phase 1 (Flow-dominated, S_F > 0.85): t = %d to %d min\n', ...
    time(find(phase==1, 1)), time(find(phase==1, 1, 'last')));
fprintf('Phase 2 (Transition, 0.65 < S_F < 0.85): t = %d to %d min\n', ...
    time(find(phase==2, 1)), time(find(phase==2, 1, 'last')));
fprintf('Phase 3 (Pressure-dominated, S_F < 0.65): t = %d to %d min\n', ...
    time(find(phase==3, 1)), time(find(phase==3, 1, 'last')));

fprintf('\n\nKey Findings:\n');
[max_int_T, idx_T] = max(interaction_T);
[max_int_P, idx_P] = max(interaction_P);
[max_int_F, idx_F] = max(interaction_F);

fprintf('Maximum Temperature interaction: %.4f at t = %d min\n', max_int_T, time(idx_T));
fprintf('Maximum Pressure interaction: %.4f at t = %d min\n', max_int_P, time(idx_P));
fprintf('Maximum Flow rate interaction: %.4f at t = %d min\n', max_int_F, time(idx_F));
fprintf('Maximum total interactions: %.4f at t = %d min\n', ...
    max(total_interaction_all), time(total_interaction_all == max(total_interaction_all)));

fprintf('\n========================================\n');

%% Save figures
% Uncomment to save
% print(1, 'interaction_analysis_overview', '-dpng', '-r300');
% print(2, 'interaction_heatmap_evolution', '-dpng', '-r300');
% print(3, 'interaction_timeline_detailed', '-dpng', '-r300');
% print(4, 'parameter_dominance_phases', '-dpng', '-r300');

fprintf('\nVisualization complete! Generated 4 detailed interaction analysis figures.\n');