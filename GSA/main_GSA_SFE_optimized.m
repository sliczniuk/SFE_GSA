close all; clear; clc;

%% Initialize parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

%% USER INPUTS

% Number of Monte Carlo samples (can use fewer with quasi-random)
N = 1e4;

% Bootstrap settings for confidence intervals
n_bootstrap = 0;  % Number of bootstrap resamples (set to 0 to disable)
% Note: Bootstrap adds computational time (~10-20% overhead)

NAMES = {'T', 'P', 'F'};
NAMES_LATEX = {'T[$^\circ C$]', 'P[bar]', 'F[kg/s]'};

input_ranges = [30   40;
                100  200;
                3.33 6.67];

% Number of input parameters
d = size(input_ranges, 1);

%% Generate Sobol quasi-random samples (better convergence than pseudo-random)
[A, B] = sobol_sampling(input_ranges, N);

%% Time points to analyze
time_points = [5, 15, 30, 45, 60, 90, 120, 150, 180, 240, 300, 360, 420, 480, 600, 720, 840, 960, 1100, 1300, 1500];
%time_points = [1200, 1400, 1600, 2000];
time_points = 300;
n_times = length(time_points);

%% Preallocate results storage

results.time_points = time_points;
results.first_order = zeros(n_times, d);
results.total_order = zeros(n_times, d);
results.second_order = NaN(n_times, d, d);  % Second-order indices (i,j pairs)
results.third_order = NaN(n_times, 1);      % Third-order index (only for d=3)
results.interaction = zeros(n_times, d);
results.output_mean = zeros(n_times, 1);
results.output_std = zeros(n_times, 1);
results.output_var = zeros(n_times, 1);
results.sum_S = zeros(n_times, 1);
results.sum_ST = zeros(n_times, 1);
results.sum_all_indices = zeros(n_times, 1);  % sum(S) + sum(S2) + S3
results.computation_time = zeros(n_times, 1);
results.time_step_used = zeros(n_times, 1);

% Confidence intervals (if bootstrap enabled)
if n_bootstrap > 0
    results.S_CI = zeros(n_times, d, 2);   % [time, param, lower/upper]
    results.ST_CI = zeros(n_times, d, 2);
    results.S3_CI = zeros(n_times, 2);     % [time, lower/upper]
    results.bootstrap_params.n_bootstrap = n_bootstrap;
    results.bootstrap_params.confidence_level = 0.95;
else
    results.S_CI = [];
    results.ST_CI = [];
    results.S3_CI = [];
    results.bootstrap_params = struct('n_bootstrap', 0);
end

results.N = N;
results.input_ranges = input_ranges;
results.NAMES = NAMES;

%% Determine unique time steps needed and build caches
% Dynamic time step selection: larger steps for longer simulations
time_steps = select_time_steps(time_points);
unique_steps = unique(time_steps);

fprintf('Building integrator caches for time steps: [%s] minutes\n', num2str(unique_steps));

caches = containers.Map('KeyType', 'double', 'ValueType', 'any');
for ts = unique_steps
    fprintf('  Building cache for time_step = %.2f min...\n', ts);
    caches(ts) = init_Extraction_Cache(ts);
end
fprintf('All caches ready.\n\n');

%% Main analysis loop
diary('myTextLog_optimized.txt');

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║         GLOBAL SENSITIVITY ANALYSIS - TIME EVOLUTION          ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Configuration:\n');
fprintf('  • Sample size (N):        %d\n', N);
fprintf('  • Number of parameters:   %d (%s)\n', d, strjoin(NAMES, ', '));
fprintf('  • Time points to analyze: %d\n', n_times);
fprintf('  • Bootstrap resamples:    %d\n', n_bootstrap);
fprintf('  • Model evaluations/time: %d\n', (2 + d + d*(d-1)/2) * N);
fprintf('\n');

total_start = tic;

for idx = 1:n_times
    time = time_points(idx);
    time_step = time_steps(idx);
    cache = caches(time_step);

    fprintf('\n╔════════════════════════════════════════════════════════════════╗\n');
    fprintf('║ TIME POINT %d/%d: t = %.0f min                                \n', idx, n_times, time);
    fprintf('╚════════════════════════════════════════════════════════════════╝\n');

    %% Create model function using cached data
    model = @(YY) Simulate_Extraction_Cached(YY, time, time_step, cache);
    model_func = @(AA) applyToEachRowOptimized(model, AA);

    %% Run Sobol analysis
    tic
    if n_bootstrap > 0
        [first_order, total_order, second_order, third_order, YA, YB, CI] = ...
            Sen_Saltelli(A, B, model_func, 'Janon-Monod', 'Jansen', 'Bootstrap', n_bootstrap);
    else
        [first_order, total_order, second_order, third_order, YA, YB, CI] = ...
            Sen_Saltelli(A, B, model_func, 'Janon-Monod', 'Jansen');
    end
    analysis_time = toc;

    %% Compute statistics
    all_outputs = [YA; YB];
    interaction = total_order - first_order;

    %% Store results
    results.first_order(idx, :)     = first_order;
    results.total_order(idx, :)     = total_order;
    results.second_order(idx, :, :) = second_order;
    results.third_order(idx)        = third_order;
    results.interaction(idx, :)     = interaction;
    results.output_mean(idx)        = mean(all_outputs);
    results.output_std(idx)         = std(all_outputs);
    results.output_var(idx)         = var(all_outputs);
    results.sum_S(idx)              = sum(first_order);
    results.sum_ST(idx)             = sum(total_order);

    % Calculate total variance accounted for
    S2_mat                       = squeeze(second_order);
    sum_S2                       = S2_mat(1,2) + S2_mat(1,3) + S2_mat(2,3);
    results.sum_all_indices(idx) = sum(first_order) + sum_S2 + third_order;

    results.computation_time(idx) = analysis_time;
    results.time_step_used(idx) = time_step;

    % Store confidence intervals if available
    if n_bootstrap > 0 && ~isempty(CI.S_CI)
        results.S_CI(idx, :, :) = CI.S_CI;
        results.ST_CI(idx, :, :) = CI.ST_CI;
        if d == 3 && ~isnan(CI.S3_CI(1))
            results.S3_CI(idx, :) = CI.S3_CI;
        end
    end

    %% Display results summary for this time point
    fprintf('\n--- RESULTS SUMMARY (t = %.0f min) ---\n', time);
    fprintf('Analysis time: %.2f seconds\n', analysis_time);
    fprintf('Output: Mean = %.4f ± %.4f\n', results.output_mean(idx), results.output_std(idx));

    fprintf('\nFirst-order indices:\n');
    for i = 1:d
        fprintf('  %s: %.4f', NAMES{i}, first_order(i));
        if n_bootstrap > 0 && ~isempty(CI.S_CI)
            fprintf(' [%.4f, %.4f]', CI.S_CI(i,1), CI.S_CI(i,2));
        end
        fprintf('\n');
    end

    fprintf('\nTotal-order indices:\n');
    for i = 1:d
        fprintf('  %s: %.4f', NAMES{i}, total_order(i));
        if n_bootstrap > 0 && ~isempty(CI.ST_CI)
            fprintf(' [%.4f, %.4f]', CI.ST_CI(i,1), CI.ST_CI(i,2));
        end
        fprintf('\n');
    end

    if d == 3
        fprintf('\nSecond-order: S2(T,P)=%.4f, S2(T,F)=%.4f, S2(P,F)=%.4f\n', ...
            second_order(1,2), second_order(1,3), second_order(2,3));
        fprintf('Third-order:  S3=%.4f\n', third_order);
        fprintf('Closure:      sum=%.6f (target: 1.0)\n', results.sum_all_indices(idx));
    end

    % Progress estimate
    if idx < n_times
        elapsed = toc(total_start);
        avg_time_per_point = elapsed / idx;
        remaining_time = (n_times - idx) * avg_time_per_point;
        fprintf('\n[Progress: %d/%d (%.1f%%) - Estimated time remaining: %.1f min]\n', ...
            idx, n_times, 100*idx/n_times, remaining_time/60);
    end

end

total_time = toc(total_start);
results.total_computation_time = total_time;

diary('off');

%% Save results
results_filename = sprintf('GSA_results_N%d_%s.mat', N, datestr(now, 'yyyy-mm-dd_HH-MM'));
%save(results_filename, 'results', 'A', 'B');

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                    ANALYSIS COMPLETE                           ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Computation Summary:\n');
fprintf('  • Total time:           %.2f minutes (%.1f hours)\n', total_time/60, total_time/3600);
fprintf('  • Time per point:       %.2f minutes\n', total_time/n_times/60);
fprintf('  • Results file:         %s\n', results_filename);
fprintf('  • Total model calls:    %d\n', (2 + d + d*(d-1)/2) * N * n_times);
fprintf('  • Average speed:        %.0f evals/sec\n', (2 + d + d*(d-1)/2) * N * n_times / total_time);
fprintf('\n');

%% Display summary table
fprintf('\n\nSUMMARY TABLE:\n');
fprintf('%-8s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s |\n', ...
    'Time', 'S_T', 'S_P', 'S_F', 'ST_T', 'ST_P', 'ST_F', 'Sum_S');
fprintf('%s\n', repmat('-', 1, 90));
for idx = 1:n_times
    fprintf('%-8.0f | %-8.4f | %-8.4f | %-8.4f | %-8.4f | %-8.4f | %-8.4f | %-8.4f | \n', ...
        results.time_points(idx), ...
        results.first_order(idx, 1), results.first_order(idx, 2), results.first_order(idx, 3), ...
        results.total_order(idx, 1), results.total_order(idx, 2), results.total_order(idx, 3), ...
        results.sum_S(idx));
end

% Display confidence intervals if available
if n_bootstrap > 0 && ~isempty(results.S_CI)
    fprintf('\n\nFIRST-ORDER INDICES WITH 95%% CONFIDENCE INTERVALS:\n');
    fprintf('%-8s | %-20s | %-20s | %-20s |\n', 'Time', 'T [CI]', 'P [CI]', 'F [CI]');
    fprintf('%s\n', repmat('-', 1, 80));
    for idx = 1:n_times
        fprintf('%-8.0f | %.4f [%.4f,%.4f] | %.4f [%.4f,%.4f] | %.4f [%.4f,%.4f] |\n', ...
            results.time_points(idx), ...
            results.first_order(idx, 1), results.S_CI(idx, 1, 1), results.S_CI(idx, 1, 2), ...
            results.first_order(idx, 2), results.S_CI(idx, 2, 1), results.S_CI(idx, 2, 2), ...
            results.first_order(idx, 3), results.S_CI(idx, 3, 1), results.S_CI(idx, 3, 2));
    end

    fprintf('\n\nTOTAL-ORDER INDICES WITH 95%% CONFIDENCE INTERVALS:\n');
    fprintf('%-8s | %-20s | %-20s | %-20s |\n', 'Time', 'T [CI]', 'P [CI]', 'F [CI]');
    fprintf('%s\n', repmat('-', 1, 80));
    for idx = 1:n_times
        fprintf('%-8.0f | %.4f [%.4f,%.4f] | %.4f [%.4f,%.4f] | %.4f [%.4f,%.4f] |\n', ...
            results.time_points(idx), ...
            results.total_order(idx, 1), results.ST_CI(idx, 1, 1), results.ST_CI(idx, 1, 2), ...
            results.total_order(idx, 2), results.ST_CI(idx, 2, 1), results.ST_CI(idx, 2, 2), ...
            results.total_order(idx, 3), results.ST_CI(idx, 3, 1), results.ST_CI(idx, 3, 2));
    end
end

fprintf('\n\nSECOND-ORDER INDICES TABLE:\n');
fprintf('%-8s | %-10s | %-10s | %-10s |\n', 'Time', 'S2(T,P)', 'S2(T,F)', 'S2(P,F)');
fprintf('%s\n', repmat('-', 1, 50));
for idx = 1:n_times
    S2_mat = squeeze(results.second_order(idx, :, :));
    fprintf('%-8.0f | %-10.4f | %-10.4f | %-10.4f |\n', ...
        results.time_points(idx), ...
        S2_mat(1,2), S2_mat(1,3), S2_mat(2,3));
end

if d == 3
    fprintf('\n\nTHIRD-ORDER INDEX & CLOSURE TABLE:\n');
    fprintf('%-8s | %-12s | %-12s |\n', 'Time', 'S3(T,P,F)', 'Total(S+S2+S3)');
    fprintf('%s\n', repmat('-', 1, 40));
    for idx = 1:n_times
        fprintf('%-8.0f | %-12.6f | %-12.6f |\n', ...
            results.time_points(idx), ...
            results.third_order(idx), ...
            results.sum_all_indices(idx));
    end
end

%% Export results to CSV
%export_to_csv(results, sprintf('GSA_results_N%d.csv', N));

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function export_to_csv(results, filename)
% EXPORT_TO_CSV Export Sobol sensitivity results to CSV file
%   Creates a comprehensive CSV file with all sensitivity indices

    fprintf('\nExporting results to CSV: %s\n', filename);

    % Create table with main results
    T = table();
    T.Time = results.time_points(:);

    % First-order indices
    for i = 1:length(results.NAMES)
        T.(sprintf('S_%s', results.NAMES{i})) = results.first_order(:, i);
    end

    % Total-order indices
    for i = 1:length(results.NAMES)
        T.(sprintf('ST_%s', results.NAMES{i})) = results.total_order(:, i);
    end

    % Interaction effects
    for i = 1:length(results.NAMES)
        T.(sprintf('Interaction_%s', results.NAMES{i})) = results.interaction(:, i);
    end

    % Second-order indices (for d=3)
    if size(results.second_order, 2) == 3
        S2_TP = squeeze(results.second_order(:, 1, 2));
        S2_TF = squeeze(results.second_order(:, 1, 3));
        S2_PF = squeeze(results.second_order(:, 2, 3));
        T.S2_TP = S2_TP;
        T.S2_TF = S2_TF;
        T.S2_PF = S2_PF;
    end

    % Third-order index
    if ~isempty(results.third_order) && ~all(isnan(results.third_order))
        T.S3_TPF = results.third_order;
    end

    % Output statistics
    T.Mean_Yield = results.output_mean;
    T.Std_Yield = results.output_std;
    T.Var_Yield = results.output_var;

    % Sums
    T.Sum_S = results.sum_S;
    T.Sum_ST = results.sum_ST;
    if ~isempty(results.sum_all_indices)
        T.Sum_All = results.sum_all_indices;
    end

    % Confidence intervals (if available)
    if ~isempty(results.S_CI)
        for i = 1:length(results.NAMES)
            T.(sprintf('S_%s_Lower', results.NAMES{i})) = results.S_CI(:, i, 1);
            T.(sprintf('S_%s_Upper', results.NAMES{i})) = results.S_CI(:, i, 2);
        end
        for i = 1:length(results.NAMES)
            T.(sprintf('ST_%s_Lower', results.NAMES{i})) = results.ST_CI(:, i, 1);
            T.(sprintf('ST_%s_Upper', results.NAMES{i})) = results.ST_CI(:, i, 2);
        end
    end

    % Computation time
    T.Computation_Time_s = results.computation_time;

    % Write to CSV
    writetable(T, filename);
    fprintf('Successfully exported %d time points to %s\n', height(T), filename);
end

function time_steps = select_time_steps(time_points)
% SELECT_TIME_STEPS Dynamically select time steps based on simulation length
%   Shorter simulations use smaller time steps for accuracy.
%   Longer simulations use larger time steps for efficiency.
%
%   Rule of thumb: aim for ~50-200 integration steps per simulation

    time_steps = zeros(size(time_points));

    for i = 1:length(time_points)
        t = time_points(i);

        if t <= 30
            % Short simulations: 0.5 min steps (~10-60 steps)
            time_steps(i) = 0.5;
        elseif t <= 150
            % Medium simulations: 1 min steps (~30-150 steps)
            time_steps(i) = 1;
        elseif t <= 300
            % Long simulations: 2 min steps (~75-300 steps)
            time_steps(i) = 2;
        elseif t <= 600
            % Long simulations: 2 min steps (~75-300 steps)
            time_steps(i) = 2;
        elseif t <= 3000
            % Very long simulations: 5 min steps (~120-240 steps)
            time_steps(i) = 5;
        else
            % Extra long simulations: 10 min steps (~150-200 steps)
            time_steps(i) = 10;
        end
    end
end

function [A, B] = sobol_sampling(input_ranges, N)
% SOBOL_SAMPLING Generate Sobol quasi-random matrices A and B
%   Uses Sobol sequences for better space-filling properties and
%   faster convergence compared to pseudo-random sampling.

    k = size(input_ranges, 1);

    % Create Sobol sequence generator
    sob = sobolset(2*k, 'Skip', 1000, 'Leap', 100);
    sob = scramble(sob, 'MatousekAffineOwen');

    % Generate 2*k columns: first k for A, next k for B
    samples = net(sob, N);

    A01 = samples(:, 1:k);
    B01 = samples(:, k+1:2*k);

    % Scale to desired ranges
    A = scale_to_ranges(A01, input_ranges);
    B = scale_to_ranges(B01, input_ranges);
end

function X = scale_to_ranges(X01, ranges)
% Scale [0,1] samples to user-defined ranges
    mins = ranges(:,1)';
    maxs = ranges(:,2)';
    X = X01 .* (maxs - mins) + mins;
end
