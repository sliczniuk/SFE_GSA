% CONVERGENCE ANALYSIS FOR SOBOL INDICES
% Identifies optimal number of samples by testing convergence
% across different sample sizes
%
% Purpose: Determine minimum N needed for stable Sobol indices
% Method: Run GSA with increasing N and track index stability

close all; clear; clc;

%% Initialize parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

%% CONFIGURATION

% Sample sizes to test (logarithmic spacing recommended)
N_values = [2500, 5000, 10000, 12000];
% For quick test, use: N_values = [100, 250, 500, 1000];

% Number of repeated runs per N (to assess stability)
n_repeats = 3;  % Increase to 5-10 for more robust estimates

% Time point to analyze (choose representative point)
analysis_time = 300;  % minutes
time_step = 5;  % integration time step

% Parameter ranges
NAMES = {'T', 'P', 'F'};
input_ranges = [30   40;      % Temperature [°C]
                100  200;      % Pressure [bar]
                3.33 6.67];    % Flow rate [kg/s]

d = size(input_ranges, 1);

% Convergence criteria
convergence_threshold = 0.01;  % 1% relative change

%% Preallocate storage

n_tests = length(N_values);
results = struct();

% Store indices for each N and repeat
results.N_values = N_values;
results.n_repeats = n_repeats;
results.first_order = zeros(n_tests, n_repeats, d);
results.total_order = zeros(n_tests, n_repeats, d);
results.second_order = NaN(n_tests, n_repeats, d, d);
results.third_order = NaN(n_tests, n_repeats);
results.computation_time = zeros(n_tests, n_repeats);
results.sum_S = zeros(n_tests, n_repeats);
results.sum_ST = zeros(n_tests, n_repeats);
results.sum_all = zeros(n_tests, n_repeats);

% Statistics across repeats
results.S_mean = zeros(n_tests, d);
results.S_std = zeros(n_tests, d);
results.ST_mean = zeros(n_tests, d);
results.ST_std = zeros(n_tests, d);
results.S3_mean = zeros(n_tests, 1);
results.S3_std = zeros(n_tests, 1);

%% Initialize cache for model evaluations
fprintf('Initializing simulation cache...\n');
cache = init_Extraction_Cache(time_step);
fprintf('Cache ready.\n\n');

%% Main convergence analysis loop

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║           SOBOL INDICES CONVERGENCE ANALYSIS                   ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Configuration:\n');
fprintf('  • Sample sizes to test: [%s]\n', num2str(N_values));
fprintf('  • Repeats per N:        %d\n', n_repeats);
fprintf('  • Analysis time:        %.0f min\n', analysis_time);
fprintf('  • Parameters:           %d (%s)\n', d, strjoin(NAMES, ', '));
fprintf('  • Convergence target:   %.2f%% relative change\n', convergence_threshold*100);
fprintf('\n');

total_start = tic;

for i = 1:n_tests
    N = N_values(i);

    fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
    fprintf('║ SAMPLE SIZE TEST %d/%d: N = %d\n', i, n_tests, N);
    fprintf('╚════════════════════════════════════════════════════════════════╝\n');

    % Run multiple times to assess stability
    for r = 1:n_repeats
        fprintf('\n[Repeat %d/%d]\n', r, n_repeats);

        % Generate new Sobol samples for each repeat
        [A, B] = sobol_sampling(input_ranges, N);

        % Create model function
        model = @(YY) Simulate_Extraction_Cached(YY, analysis_time, time_step, cache);
        model_func = @(AA) applyToEachRowOptimized(model, AA);

        % Run Sobol analysis (no bootstrap)
        tic_analysis = tic;
        [first_order, total_order, second_order, third_order, ~, ~, ~] = ...
            Sen_Saltelli(A, B, model_func, 'Janon-Monod', 'Jansen');
        comp_time = toc(tic_analysis);

        % Store results
        results.first_order(i, r, :) = first_order;
        results.total_order(i, r, :) = total_order;
        results.second_order(i, r, :, :) = second_order;
        results.third_order(i, r) = third_order;
        results.computation_time(i, r) = comp_time;
        results.sum_S(i, r) = sum(first_order);
        results.sum_ST(i, r) = sum(total_order);

        % Calculate total variance accounted
        S2_mat = squeeze(second_order);
        sum_S2 = S2_mat(1,2) + S2_mat(1,3) + S2_mat(2,3);
        results.sum_all(i, r) = sum(first_order) + sum_S2 + third_order;

        % Quick display
        fprintf('  Time: %.2f s | sum(S)=%.4f | S3=%.6f | Closure=%.6f\n', ...
            comp_time, sum(first_order), third_order, results.sum_all(i, r));
    end

    % Compute statistics across repeats
    results.S_mean(i, :) = mean(squeeze(results.first_order(i, :, :)), 1);
    results.S_std(i, :) = std(squeeze(results.first_order(i, :, :)), 0, 1);
    results.ST_mean(i, :) = mean(squeeze(results.total_order(i, :, :)), 1);
    results.ST_std(i, :) = std(squeeze(results.total_order(i, :, :)), 0, 1);
    results.S3_mean(i) = mean(results.third_order(i, :));
    results.S3_std(i) = std(results.third_order(i, :));

    % Display statistics for this N
    fprintf('\n--- Statistics for N = %d (across %d repeats) ---\n', N, n_repeats);
    fprintf('First-order indices (mean ± std):\n');
    for p = 1:d
        fprintf('  %s: %.4f ± %.4f (CV=%.1f%%)\n', NAMES{p}, ...
            results.S_mean(i, p), results.S_std(i, p), ...
            100*results.S_std(i, p)/results.S_mean(i, p));
    end
    fprintf('Third-order: %.6f ± %.6f\n', results.S3_mean(i), results.S3_std(i));
    fprintf('Avg computation time: %.2f ± %.2f seconds\n', ...
        mean(results.computation_time(i, :)), std(results.computation_time(i, :)));

    % Check convergence (if not first iteration)
    if i > 1
        fprintf('\n--- Convergence check (relative change from N=%d to N=%d) ---\n', ...
            N_values(i-1), N);

        converged_all = true;
        for p = 1:d
            rel_change_S = abs(results.S_mean(i, p) - results.S_mean(i-1, p)) / ...
                          (results.S_mean(i-1, p) + 1e-10);
            rel_change_ST = abs(results.ST_mean(i, p) - results.ST_mean(i-1, p)) / ...
                           (results.ST_mean(i-1, p) + 1e-10);

            fprintf('  %s: ΔS=%.2f%%, ΔST=%.2f%%', NAMES{p}, ...
                rel_change_S*100, rel_change_ST*100);

            if rel_change_S < convergence_threshold && rel_change_ST < convergence_threshold
                fprintf(' [CONVERGED]\n');
            else
                fprintf(' [NOT CONVERGED]\n');
                converged_all = false;
            end
        end

        % Third-order convergence
        rel_change_S3 = abs(results.S3_mean(i) - results.S3_mean(i-1)) / ...
                       (abs(results.S3_mean(i-1)) + 1e-10);
        fprintf('  S3: Δ=%.2f%%', rel_change_S3*100);
        if rel_change_S3 < convergence_threshold
            fprintf(' [CONVERGED]\n');
        else
            fprintf(' [NOT CONVERGED]\n');
            converged_all = false;
        end

        if converged_all
            fprintf('\n*** ALL INDICES CONVERGED AT N = %d ***\n', N);
        end
    end

    % Progress estimate
    if i < n_tests
        elapsed = toc(total_start);
        avg_time_per_test = elapsed / i;
        remaining_time = (n_tests - i) * avg_time_per_test;
        fprintf('\n[Overall Progress: %d/%d (%.1f%%) - ETA: %.1f min]\n', ...
            i, n_tests, 100*i/n_tests, remaining_time/60);
    end
end

total_time = toc(total_start);

%% Analysis and Recommendations

fprintf('\n\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                    CONVERGENCE ANALYSIS COMPLETE               ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

% Find optimal N
fprintf('=== OPTIMAL SAMPLE SIZE RECOMMENDATION ===\n\n');

% Identify where all indices converged
optimal_N_idx = NaN;
for i = 2:n_tests
    converged = true;

    % Check all first-order indices
    for p = 1:d
        rel_change_S = abs(results.S_mean(i, p) - results.S_mean(i-1, p)) / ...
                      (results.S_mean(i-1, p) + 1e-10);
        rel_change_ST = abs(results.ST_mean(i, p) - results.ST_mean(i-1, p)) / ...
                       (results.ST_mean(i-1, p) + 1e-10);

        if rel_change_S >= convergence_threshold || rel_change_ST >= convergence_threshold
            converged = false;
            break;
        end
    end

    % Check third-order
    rel_change_S3 = abs(results.S3_mean(i) - results.S3_mean(i-1)) / ...
                   (abs(results.S3_mean(i-1)) + 1e-10);
    if rel_change_S3 >= convergence_threshold
        converged = false;
    end

    if converged
        optimal_N_idx = i;
        break;
    end
end

if ~isnan(optimal_N_idx)
    fprintf('Recommended N: %d\n', N_values(optimal_N_idx));
    fprintf('  • All indices converged within %.2f%% tolerance\n', convergence_threshold*100);
    fprintf('  • Avg computation time: %.2f seconds\n', mean(results.computation_time(optimal_N_idx, :)));
    fprintf('  • Required model evaluations: %d\n', (2 + d + d*(d-1)/2) * N_values(optimal_N_idx));
else
    fprintf('WARNING: Convergence not achieved within tested range.\n');
    fprintf('  • Recommend testing larger N values: [%d, %d, %d]\n', ...
        2*N_values(end), 5*N_values(end), 10*N_values(end));
    fprintf('  • Current best N: %d (lowest variability)\n', N_values(end));
end

fprintf('\n');
fprintf('Total analysis time: %.2f minutes\n', total_time/60);

%% Summary Table

fprintf('\n\n=== SUMMARY TABLE ===\n\n');
fprintf('%-8s | ', 'N');
for p = 1:d
    fprintf('%-12s | ', sprintf('S_%s (CV%%)', NAMES{p}));
end
fprintf('%-12s | %-12s | %-10s |\n', 'S3 (mean)', 'Closure', 'Time (s)');
fprintf('%s\n', repmat('-', 1, 100));

for i = 1:n_tests
    fprintf('%-8d | ', N_values(i));

    % First-order with CV
    for p = 1:d
        cv = 100 * results.S_std(i, p) / results.S_mean(i, p);
        fprintf('%.4f(%.1f%%) | ', results.S_mean(i, p), cv);
    end

    % Third-order and closure
    fprintf('%.6f     | ', results.S3_mean(i));
    fprintf('%.6f     | ', mean(results.sum_all(i, :)));
    fprintf('%.2f       |\n', mean(results.computation_time(i, :)));
end

%% Visualization

figure('Position', [100, 100, 1400, 900]);

% Plot 1: First-order index convergence
subplot(2, 3, 1);
for p = 1:d
    errorbar(N_values, results.S_mean(:, p), results.S_std(:, p), ...
        'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', NAMES{p});
    hold on;
end
set(gca, 'XScale', 'log');
xlabel('Sample Size (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('First-order Index', 'FontSize', 12, 'FontWeight', 'bold');
title('First-order Convergence', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% Plot 2: Total-order index convergence
subplot(2, 3, 2);
for p = 1:d
    errorbar(N_values, results.ST_mean(:, p), results.ST_std(:, p), ...
        's-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', NAMES{p});
    hold on;
end
set(gca, 'XScale', 'log');
xlabel('Sample Size (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Total-order Index', 'FontSize', 12, 'FontWeight', 'bold');
title('Total-order Convergence', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% Plot 3: Third-order convergence
subplot(2, 3, 3);
errorbar(N_values, results.S3_mean, results.S3_std, ...
    'd-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.5 0 0.5]);
set(gca, 'XScale', 'log');
xlabel('Sample Size (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Third-order Index (S3)', 'FontSize', 12, 'FontWeight', 'bold');
title('Third-order Convergence', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Plot 4: Coefficient of variation (stability metric)
subplot(2, 3, 4);
for p = 1:d
    cv = 100 * results.S_std(:, p) ./ results.S_mean(:, p);
    plot(N_values, cv, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', NAMES{p});
    hold on;
end
yline(5, 'k--', 'LineWidth', 1.5, 'DisplayName', '5% threshold');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('Sample Size (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Coefficient of Variation (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('Index Stability (CV)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 11);

% Plot 5: Closure property check
subplot(2, 3, 5);
closure_mean = mean(results.sum_all, 2);
closure_std = std(results.sum_all, 0, 2);
errorbar(N_values, closure_mean, closure_std, ...
    '^-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0.5 0]);
hold on;
yline(1.0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Perfect closure');
fill([N_values(1) N_values(end) N_values(end) N_values(1)], ...
    [0.99 0.99 1.01 1.01], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
set(gca, 'XScale', 'log');
xlabel('Sample Size (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('sum(S) + sum(S2) + S3', 'FontSize', 12, 'FontWeight', 'bold');
title('Closure Property Check', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0.95 1.05]);
grid on;
set(gca, 'FontSize', 11);

% Plot 6: Computation time scaling
subplot(2, 3, 6);
time_mean = mean(results.computation_time, 2);
time_std = std(results.computation_time, 0, 2);
errorbar(N_values, time_mean, time_std, ...
    'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.8 0.4 0]);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('Sample Size (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Computation Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Computational Cost', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Overall title
sgtitle('Sobol Indices Convergence Analysis', 'FontSize', 16, 'FontWeight', 'bold');

%% Save results
save_filename = sprintf('convergence_results_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM'));
save(save_filename, 'results');

fprintf('\n\nResults saved to: %s\n', save_filename);
fprintf('Figure displayed - close to continue or export manually.\n');

%% Helper function
function [A, B] = sobol_sampling(input_ranges, N)
    % Generate Sobol sequences for matrices A and B
    k = size(input_ranges, 1);

    % Create Sobol sequence
    p = sobolset(2*k, 'Skip', 1e3, 'Leap', 1e2);
    p = scramble(p, 'MatousekAffineOwen');

    % Generate samples
    samples = net(p, N);

    A01 = samples(:, 1:k);
    B01 = samples(:, k+1:2*k);

    % Scale to ranges
    A = scale_to_ranges(A01, input_ranges);
    B = scale_to_ranges(B01, input_ranges);
end

function X = scale_to_ranges(X01, ranges)
    % Scale [0,1] samples to user-defined ranges
    mins = ranges(:,1)';
    maxs = ranges(:,2)';
    X = X01 .* (maxs - mins) + mins;
end
