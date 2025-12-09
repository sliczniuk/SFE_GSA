% Test script to compare different Sobol index estimators
% Compares: Saltelli, Janon-Monod, Sobol (S) and Jansen, Sobol, Saltelli (ST)

close all; clear; clc;

%% Initialize
if isempty(gcp('nocreate'))
    parpool;
end

%% Test parameters
N = 5000;  % Number of samples (use smaller N for quick test)

NAMES = {'T', 'P', 'F'};

input_ranges = [30   40;
                100  200;
                3.33 6.67];

d = size(input_ranges, 1);

%% Generate Sobol samples
[A, B] = sobol_sampling(input_ranges, N);

%% Select time point for comparison
time = 300;  % 300 minutes - medium length simulation
time_step = 2;  % 2 min steps

fprintf('========================================\n');
fprintf('Sobol Estimator Comparison Test\n');
fprintf('Time = %.0f min, N = %d samples\n', time, N);
fprintf('========================================\n\n');

%% Initialize cache
fprintf('Initializing extraction model cache...\n');
cache = init_Extraction_Cache(time_step);

%% Create model function
model = @(YY) Simulate_Extraction_Cached(YY, time, time_step, cache);
model_func = @(AA) applyToEachRowOptimized(model, AA);

%% Run comparison with all estimators
fprintf('\nRunning Sobol analysis with all estimators...\n\n');
[S, ST, estimators] = Sen_Sobol_Estimators(A, B, model_func);

%% Display detailed comparison
fprintf('\n========================================\n');
fprintf('DETAILED RESULTS\n');
fprintf('========================================\n');

fprintf('\n--- First-Order Indices (S) ---\n');
fprintf('%-12s | %-12s | %-12s | %-12s\n', 'Parameter', 'Saltelli', 'Janon-Monod', 'Sobol');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:d
    fprintf('%-12s | %12.6f | %12.6f | %12.6f\n', ...
        NAMES{i}, estimators.S_Saltelli(i), estimators.S_Janon(i), estimators.S_Sobol(i));
end
fprintf('%-12s | %12.6f | %12.6f | %12.6f\n', ...
    'Sum', sum(estimators.S_Saltelli), sum(estimators.S_Janon), sum(estimators.S_Sobol));

fprintf('\n--- Total-Effect Indices (ST) ---\n');
fprintf('%-12s | %-12s | %-12s | %-12s\n', 'Parameter', 'Jansen', 'Sobol', 'Saltelli');
fprintf('%s\n', repmat('-', 1, 55));
for i = 1:d
    fprintf('%-12s | %12.6f | %12.6f | %12.6f\n', ...
        NAMES{i}, estimators.ST_Jansen(i), estimators.ST_Sobol(i), estimators.ST_Saltelli(i));
end
fprintf('%-12s | %12.6f | %12.6f | %12.6f\n', ...
    'Sum', sum(estimators.ST_Jansen), sum(estimators.ST_Sobol), sum(estimators.ST_Saltelli));

fprintf('\n--- Interaction Terms (ST - S) ---\n');
fprintf('Using Jansen (ST) and Janon-Monod (S):\n');
interaction = estimators.ST_Jansen - estimators.S_Janon;
for i = 1:d
    status = '';
    if interaction(i) < 0
        status = ' (WARNING: negative!)';
    end
    fprintf('  %s: %.6f%s\n', NAMES{i}, interaction(i), status);
end

fprintf('\n--- Estimator Agreement (std across methods) ---\n');
S_all = [estimators.S_Saltelli; estimators.S_Janon; estimators.S_Sobol];
ST_all = [estimators.ST_Jansen; estimators.ST_Sobol; estimators.ST_Saltelli];
fprintf('First-order std:   [%s]\n', sprintf('%.6f ', std(S_all)));
fprintf('Total-effect std:  [%s]\n', sprintf('%.6f ', std(ST_all)));

fprintf('\n--- Output Statistics ---\n');
fprintf('Mean output (f0): %.6f\n', estimators.f0);
fprintf('Variance: %.6f\n', estimators.VarY);

%% Save results
results.estimators = estimators;
results.S = S;
results.ST = ST;
results.time = time;
results.N = N;
results.NAMES = NAMES;

save(sprintf('estimator_comparison_t%d_N%d.mat', time, N), 'results');
fprintf('\nResults saved to estimator_comparison_t%d_N%d.mat\n', time, N);

%% Helper function
function [A, B] = sobol_sampling(input_ranges, N)
    k = size(input_ranges, 1);
    sob = sobolset(2*k, 'Skip', 1000, 'Leap', 100);
    sob = scramble(sob, 'MatousekAffineOwen');
    samples = net(sob, N);
    A01 = samples(:, 1:k);
    B01 = samples(:, k+1:2*k);
    mins = input_ranges(:,1)';
    maxs = input_ranges(:,2)';
    A = A01 .* (maxs - mins) + mins;
    B = B01 .* (maxs - mins) + mins;
end
