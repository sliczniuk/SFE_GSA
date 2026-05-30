close all; clear; clc;

%% Resolve paths without modifying the original project files
this_file = mfilename('fullpath');
this_folder = fileparts(this_file);
parent_folder = fileparts(this_folder);

addpath(this_folder);
addpath(parent_folder);

original_folder = pwd;
cleanup_obj = onCleanup(@() cd(original_folder));
cd(parent_folder);

%% Initialize parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

%% USER INPUTS
N = 1e4;
n_bootstrap = 0;

NAMES = {'T', 'P', 'F'};
NAMES_LATEX = {'T[$^\circ C$]', 'P[bar]', 'F[kg/s]'};

input_ranges = [30   40;
                100  200;
                3.33 6.67];

d = size(input_ranges, 1);

%% Generate Sobol quasi-random samples
[A, B] = sobol_sampling(input_ranges, N);

%% Time points to analyze
%time_points = [5, 15, 30, 45, 60, 90, 120, 150, 180, 240, 300, ...
%               360, 420, 480, 600, 720, 840, 960, 1100, 1300, 1500];
time_points = [15, 60];
n_times = numel(time_points);

%% Use one integration grid for the full trajectory
candidate_time_steps = select_time_steps(time_points);
time_step = min(candidate_time_steps);
time_step = 1;

step_idx = round(time_points ./ time_step);
tolerance = max(1e-10, 1e-9 * max(1, max(time_points)));
if any(abs(step_idx .* time_step - time_points) > tolerance)
    error('time_points must align with the selected time_step %.12g.', time_step);
end

%% Preallocate results storage
results.time_points = time_points;
results.first_order = zeros(n_times, d);
results.total_order = zeros(n_times, d);
results.second_order = NaN(n_times, d, d);
results.third_order = NaN(n_times, 1);
results.interaction = zeros(n_times, d);
results.output_mean = zeros(n_times, 1);
results.output_std = zeros(n_times, 1);
results.output_var = zeros(n_times, 1);
results.sum_S = zeros(n_times, 1);
results.sum_ST = zeros(n_times, 1);
results.sum_all_indices = zeros(n_times, 1);
results.computation_time = zeros(n_times, 1);
results.time_step_used = time_step * ones(n_times, 1);
results.single_pass = true;
results.N = N;
results.input_ranges = input_ranges;
results.NAMES = NAMES;
results.NAMES_LATEX = NAMES_LATEX;

if n_bootstrap > 0
    results.S_CI = NaN(n_times, d, 2);
    results.ST_CI = NaN(n_times, d, 2);
    results.S2_CI = NaN(n_times, d, d, 2);
    results.S3_CI = NaN(n_times, 2);
    results.bootstrap_params.n_bootstrap = n_bootstrap;
    results.bootstrap_params.confidence_level = 0.95;
else
    results.S_CI = [];
    results.ST_CI = [];
    results.S2_CI = [];
    results.S3_CI = [];
    results.bootstrap_params = struct('n_bootstrap', 0);
end

%% Build one cache and run one time-series Sobol analysis
fprintf('Building integrator cache for single time_step = %.2f min\n', time_step);
cache = init_Extraction_Cache(time_step);

diary_file = fullfile(this_folder, 'myTextLog_single_pass.txt');
diary(diary_file);

fprintf('\n');
fprintf('================================================================\n');
fprintf(' GLOBAL SENSITIVITY ANALYSIS - SINGLE-PASS TIME EVOLUTION\n');
fprintf('================================================================\n\n');
fprintf('Configuration:\n');
fprintf('  Sample size (N):        %d\n', N);
fprintf('  Number of parameters:   %d (%s)\n', d, strjoin(NAMES, ', '));
fprintf('  Time points to analyze: %d\n', n_times);
fprintf('  Single time step:       %.4g min\n', time_step);
fprintf('  Bootstrap resamples:    %d\n', n_bootstrap);

old_evals = (2 + d + d * (d - 1) / 2) * N * n_times;
new_evals = (2 + d + d * (d - 1) / 2) * N;
fprintf('  Old model row evals:    %d\n', old_evals);
fprintf('  New model row evals:    %d\n', new_evals);
fprintf('  Evaluation reduction:   %.1fx\n\n', old_evals / new_evals);

model = @(YY) Simulate_Extraction_Cached_Trajectory(YY, time_points, time_step, cache);
model_func = @(AA) applyToEachRowMatrix(model, AA);

total_start = tic;
if n_bootstrap > 0
    [first_order, total_order, second_order, third_order, YA, YB, CI, evals] = ...
        Sen_Saltelli_TimeSeries(A, B, model_func, 'Janon-Monod', 'Jansen', 'Bootstrap', n_bootstrap);
else
    [first_order, total_order, second_order, third_order, YA, YB, CI, evals] = ...
        Sen_Saltelli_TimeSeries(A, B, model_func, 'Janon-Monod', 'Jansen');
end
analysis_time = toc(total_start);

%% Store time-indexed results
all_outputs = [YA; YB];
interaction = total_order - first_order;

results.first_order = first_order;
results.total_order = total_order;
results.second_order = second_order;
results.third_order = third_order;
results.interaction = interaction;
results.output_mean = mean(all_outputs, 1)';
results.output_std = std(all_outputs, 0, 1)';
results.output_var = var(all_outputs, 0, 1)';
results.sum_S = sum(first_order, 2);
results.sum_ST = sum(total_order, 2);
results.computation_time(:) = analysis_time / n_times;
results.single_pass_analysis_time = analysis_time;
results.model_evaluations_old = old_evals;
results.model_evaluations_new = evals.total_model_rows;
results.evaluation_reduction = old_evals / evals.total_model_rows;

if d == 3
    results.sum_all_indices = sum(first_order, 2) + ...
        squeeze(second_order(:, 1, 2)) + ...
        squeeze(second_order(:, 1, 3)) + ...
        squeeze(second_order(:, 2, 3)) + third_order;
end

if n_bootstrap > 0
    results.S_CI = CI.S_CI;
    results.ST_CI = CI.ST_CI;
    results.S2_CI = CI.S2_CI;
    results.S3_CI = CI.S3_CI;
end

%% Validate dimensions
assert(isequal(size(results.first_order), [n_times, d]), 'Unexpected first_order dimensions.');
assert(isequal(size(results.total_order), [n_times, d]), 'Unexpected total_order dimensions.');
assert(isequal(size(results.second_order), [n_times, d, d]), 'Unexpected second_order dimensions.');
assert(isequal(size(YA), [N, n_times]), 'Unexpected YA dimensions.');
assert(isequal(size(YB), [N, n_times]), 'Unexpected YB dimensions.');

diary('off');

%% Save results in the single-pass folder
results_filename = fullfile(this_folder, sprintf('GSA_single_pass_results_N%d_%s.mat', N, datestr(now, 'yyyy-mm-dd_HH-MM')));
save(results_filename, 'results', 'A', 'B', 'YA', 'YB');

fprintf('\n');
fprintf('================================================================\n');
fprintf(' ANALYSIS COMPLETE\n');
fprintf('================================================================\n\n');
fprintf('Computation Summary:\n');
fprintf('  Total time:            %.2f minutes (%.1f hours)\n', analysis_time / 60, analysis_time / 3600);
fprintf('  Results file:          %s\n', results_filename);
fprintf('  Model row evals:       %d\n', evals.total_model_rows);
fprintf('  Evaluation reduction:  %.1fx\n', results.evaluation_reduction);

fprintf('\n\nSUMMARY TABLE:\n');
fprintf('%-8s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s |\n', ...
    'Time', 'S_T', 'S_P', 'S_F', 'ST_T', 'ST_P', 'ST_F', 'Sum_S');
fprintf('%s\n', repmat('-', 1, 90));
for idx = 1:n_times
    fprintf('%-8.0f | %-8.4f | %-8.4f | %-8.4f | %-8.4f | %-8.4f | %-8.4f | %-8.4f |\n', ...
        results.time_points(idx), ...
        results.first_order(idx, 1), results.first_order(idx, 2), results.first_order(idx, 3), ...
        results.total_order(idx, 1), results.total_order(idx, 2), results.total_order(idx, 3), ...
        results.sum_S(idx));
end

if d == 3
    fprintf('\n\nSECOND-ORDER AND THIRD-ORDER INDICES:\n');
    fprintf('%-8s | %-10s | %-10s | %-10s | %-10s | %-12s |\n', ...
        'Time', 'S2(T,P)', 'S2(T,F)', 'S2(P,F)', 'S3', 'SumAll');
    fprintf('%s\n', repmat('-', 1, 75));
    for idx = 1:n_times
        fprintf('%-8.0f | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-12.6f |\n', ...
            results.time_points(idx), ...
            results.second_order(idx, 1, 2), ...
            results.second_order(idx, 1, 3), ...
            results.second_order(idx, 2, 3), ...
            results.third_order(idx), ...
            results.sum_all_indices(idx));
    end
end

%% Helper functions
function time_steps = select_time_steps(time_points)
    time_steps = zeros(size(time_points));

    for i = 1:length(time_points)
        t = time_points(i);

        if t <= 30
            time_steps(i) = 0.5;
        elseif t <= 150
            time_steps(i) = 1;
        elseif t <= 600
            time_steps(i) = 2;
        elseif t <= 3000
            time_steps(i) = 5;
        else
            time_steps(i) = 10;
        end
    end
end

function [A, B] = sobol_sampling(input_ranges, N)
    k = size(input_ranges, 1);
    sob = sobolset(2 * k, 'Skip', 1000, 'Leap', 100);
    sob = scramble(sob, 'MatousekAffineOwen');
    samples = net(sob, N);

    A01 = samples(:, 1:k);
    B01 = samples(:, k + 1:2 * k);

    A = scale_to_ranges(A01, input_ranges);
    B = scale_to_ranges(B01, input_ranges);
end

function X = scale_to_ranges(X01, ranges)
    mins = ranges(:, 1)';
    maxs = ranges(:, 2)';
    X = X01 .* (maxs - mins) + mins;
end
