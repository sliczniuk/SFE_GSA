%% diagnostic_GSA_before_full_run.m
% Diagnostic checks before rerunning full Sobol analysis on SFE model.
%
% Purpose:
%   1. Verify Sobol estimator on benchmark functions.
%   2. Check convergence with increasing N.
%   3. Check repeated scrambled Sobol sequences.
%   4. Check model determinism.
%   5. Optionally check time-step sensitivity for the SFE model.
%
% Requirements:
%   - Sen_Saltelli.m available in path
%   - SFE model functions available if running model diagnostics:
%       init_Extraction_Cache
%       Simulate_Extraction_Cached
%       applyToEachRowOptimized

close all; clear; clc;

%% ------------------------------------------------------------------------
% USER SETTINGS
% -------------------------------------------------------------------------

run_benchmark_tests = true;
run_sfe_tests       = true;

% Representative sample sizes for convergence study
N_list = 2.^[13];  % Increase to 2^16 if affordable

% Number of independent scrambled Sobol repetitions
n_repeats = 6;

% Inputs: T, P, F
NAMES = {'T', 'P', 'F'};

input_ranges = [30   40;
                100  200;
                3.33 6.67];

d = size(input_ranges, 1);

% Representative time points for SFE diagnostics
test_time_points = [100];

% Optional manually imposed time steps for time-step sensitivity
% Leave empty to use select_time_steps().
manual_time_steps_to_test = [15, 10, 5];  % Useful for t = 1500 min

% Estimators used in your main analysis
first_order_estimator = 'Janon-Monod';
total_order_estimator = 'Jansen';

fprintf('\n');
fprintf('============================================================\n');
fprintf('       GSA / SOBOL DIAGNOSTIC SCRIPT\n');
fprintf('============================================================\n\n');

fprintf('Inputs: %s\n', strjoin(NAMES, ', '));
fprintf('Dimension d = %d\n', d);
fprintf('Sample sizes tested: %s\n', num2str(N_list));
fprintf('Repeated scramblings: %d\n\n', n_repeats);

%% ------------------------------------------------------------------------
% 1. BENCHMARK TESTS
% -------------------------------------------------------------------------

if run_benchmark_tests

    fprintf('\n============================================================\n');
    fprintf('1. BENCHMARK TESTS\n');
    fprintf('============================================================\n');

    benchmark_models = {};
    benchmark_names = {};

    % Additive model: no interactions
    benchmark_models{end+1} = @(X) X(:,1) + X(:,2) + X(:,3);
    benchmark_names{end+1} = 'Additive: y = x1 + x2 + x3';

    % Pairwise interaction model
    benchmark_models{end+1} = @(X) X(:,1) + X(:,2) + X(:,3) + X(:,1).*X(:,2);
    benchmark_names{end+1} = 'Pairwise interaction: y = x1 + x2 + x3 + x1*x2';

    % Pure multiplicative model
    benchmark_models{end+1} = @(X) X(:,1).*X(:,2).*X(:,3);
    benchmark_names{end+1} = 'Multiplicative: y = x1*x2*x3';

    % Use normalized [0,1] ranges for benchmark tests
    benchmark_ranges = repmat([0 1], d, 1);

    for m = 1:numel(benchmark_models)

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Benchmark model %d: %s\n', m, benchmark_names{m});
        fprintf('------------------------------------------------------------\n');

        Tbench = table();

        for n_idx = 1:numel(N_list)

            N = N_list(n_idx);

            [A, B] = sobol_sampling_diagnostic(benchmark_ranges, N, d, 1000 + n_idx);

            model_func = benchmark_models{m};

            [S, ST, S2, S3] = Sen_Saltelli( ...
                A, B, model_func, first_order_estimator, total_order_estimator);

            violation = max(S - ST);

            S2_pairs = get_second_order_pairs(S2);

            sum_all = sum(S) + sum(S2_pairs) + S3;

            Tbench.N(n_idx,1) = N;
            Tbench.S1(n_idx,1) = S(1);
            Tbench.S2(n_idx,1) = S(2);
            Tbench.S3_input(n_idx,1) = S(3);
            Tbench.ST1(n_idx,1) = ST(1);
            Tbench.ST2(n_idx,1) = ST(2);
            Tbench.ST3_input(n_idx,1) = ST(3);
            Tbench.S12(n_idx,1) = S2(1,2);
            Tbench.S13(n_idx,1) = S2(1,3);
            Tbench.S23(n_idx,1) = S2(2,3);
            Tbench.S123(n_idx,1) = S3;
            Tbench.Max_S_minus_ST(n_idx,1) = violation;
            Tbench.Closure(n_idx,1) = sum_all;

        end

        disp(Tbench);

        fprintf('\nInterpretation:\n');
        fprintf('  - Additive model should have S_i approximately equal to ST_i.\n');
        fprintf('  - Additive model should have second- and third-order terms near zero.\n');
        fprintf('  - Max_S_minus_ST should be close to zero or negative.\n');

    end
end

%% ------------------------------------------------------------------------
% 2. SFE MODEL DETERMINISM TEST
% -------------------------------------------------------------------------

if run_sfe_tests

    fprintf('\n============================================================\n');
    fprintf('2. SFE MODEL DETERMINISM TEST\n');
    fprintf('============================================================\n');

    % Choose a representative input point
    x_test = mean(input_ranges, 2)';

    fprintf('Testing repeated model calls at input:\n');
    for i = 1:d
        fprintf('  %s = %.6g\n', NAMES{i}, x_test(i));
    end

    for t_idx = 1:numel(test_time_points)

        time = test_time_points(t_idx);
        time_step = select_time_steps_diagnostic(time);

        fprintf('\nTime = %.0f min, time_step = %.3g min\n', time, time_step);

        cache = init_Extraction_Cache(time_step);
        model = @(YY) Simulate_Extraction_Cached(YY, time, time_step, cache);

        n_repeat_model = 5;
        y = zeros(n_repeat_model,1);

        for r = 1:n_repeat_model
            y(r) = model(x_test);
        end

        disp(table((1:n_repeat_model)', y, ...
            'VariableNames', {'Repeat', 'ModelOutput'}));

        fprintf('Max absolute difference = %.4e\n', max(abs(y - y(1))));

        if max(abs(y - y(1))) > 1e-8
            warning('Model output is not exactly deterministic at this point.');
        end
    end
end

%% ------------------------------------------------------------------------
% 3. SFE CONVERGENCE WITH N AT SELECTED TIME POINTS
% -------------------------------------------------------------------------

if run_sfe_tests

    fprintf('\n============================================================\n');
    fprintf('3. SFE CONVERGENCE WITH SAMPLE SIZE N\n');
    fprintf('============================================================\n');

    all_convergence_tables = struct();

    for t_idx = 1:numel(test_time_points)

        time = test_time_points(t_idx);
        time_step = select_time_steps_diagnostic(time);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('SFE convergence test: time = %.0f min, time_step = %.3g min\n', ...
            time, time_step);
        fprintf('------------------------------------------------------------\n');

        cache = init_Extraction_Cache(time_step);
        model = @(YY) Simulate_Extraction_Cached(YY, time, time_step, cache);
        model_func = @(AA) applyToEachRowOptimized(model, AA);

        Tconv = table();

        for n_idx = 1:numel(N_list)

            N = N_list(n_idx);

            fprintf('\nRunning N = %d...\n', N);

            [A, B] = sobol_sampling_diagnostic(input_ranges, N, d, 2000 + n_idx);

            tic;
            [S, ST, S2, S123, YA, YB] = Sen_Saltelli( ...
                A, B, model_func, first_order_estimator, total_order_estimator);
            elapsed = toc;

            S2_pairs = get_second_order_pairs(S2);
            closure = sum(S) + sum(S2_pairs) + S123;

            Tconv.N(n_idx,1) = N;

            for i = 1:d
                Tconv.(sprintf('S_%s', NAMES{i}))(n_idx,1) = S(i);
                Tconv.(sprintf('ST_%s', NAMES{i}))(n_idx,1) = ST(i);
            end

            Tconv.S2_TP(n_idx,1) = S2(1,2);
            Tconv.S2_TF(n_idx,1) = S2(1,3);
            Tconv.S2_PF(n_idx,1) = S2(2,3);
            Tconv.S123(n_idx,1) = S123;

            Tconv.Sum_S(n_idx,1) = sum(S);
            Tconv.Closure(n_idx,1) = closure;
            Tconv.Max_S_minus_ST(n_idx,1) = max(S - ST);
            Tconv.OutputMean(n_idx,1) = mean([YA; YB]);
            Tconv.OutputStd(n_idx,1) = std([YA; YB]);
            Tconv.Elapsed_s(n_idx,1) = elapsed;

            fprintf('  S  = [% .4f, % .4f, % .4f]\n', S);
            fprintf('  ST = [% .4f, % .4f, % .4f]\n', ST);
            fprintf('  max(S-ST) = %.4e\n', max(S - ST));
            fprintf('  closure   = %.6f\n', closure);

        end

        disp(Tconv);

        field_name = sprintf('t_%g', time);
        field_name = strrep(field_name, '.', '_');
        all_convergence_tables.(field_name) = Tconv;

        filename = sprintf('diagnostic_convergence_t_%g.csv', time);
        writetable(Tconv, filename);
        fprintf('Saved convergence table to %s\n', filename);

    end

    save('diagnostic_convergence_results.mat', 'all_convergence_tables');

end

%% ------------------------------------------------------------------------
% 4. REPEATED SCRAMBLED SOBOL SEQUENCES
% -------------------------------------------------------------------------

if run_sfe_tests

    fprintf('\n============================================================\n');
    fprintf('4. REPEATED SCRAMBLED SOBOL SEQUENCE TEST\n');
    fprintf('============================================================\n');

    % Use one moderate/large N for repeated scrambling
    N_repeat_test = N_list(end);

    repeated_results = struct();

    for t_idx = 1:numel(test_time_points)

        time = test_time_points(t_idx);
        time_step = select_time_steps_diagnostic(time);

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Repeated scrambling test: time = %.0f min, N = %d\n', ...
            time, N_repeat_test);
        fprintf('------------------------------------------------------------\n');

        cache = init_Extraction_Cache(time_step);
        model = @(YY) Simulate_Extraction_Cached(YY, time, time_step, cache);
        model_func = @(AA) applyToEachRowOptimized(model, AA);

        S_all = zeros(n_repeats, d);
        ST_all = zeros(n_repeats, d);
        S2_all = zeros(n_repeats, 3);
        S123_all = zeros(n_repeats, 1);
        closure_all = zeros(n_repeats, 1);
        violation_all = zeros(n_repeats, 1);

        for r = 1:n_repeats

            fprintf('  Scrambling repeat %d/%d...\n', r, n_repeats);

            [A, B] = sobol_sampling_diagnostic(input_ranges, N_repeat_test, d, 5000 + r);

            [S, ST, S2, S123] = Sen_Saltelli( ...
                A, B, model_func, first_order_estimator, total_order_estimator);

            S_all(r,:) = S(:)';
            ST_all(r,:) = ST(:)';
            S2_all(r,:) = [S2(1,2), S2(1,3), S2(2,3)];
            S123_all(r) = S123;
            closure_all(r) = sum(S) + sum(S2_all(r,:)) + S123;
            violation_all(r) = max(S - ST);

        end

        Trep = table();
        Trep.Index = {'S_T'; 'S_P'; 'S_F'; ...
                      'ST_T'; 'ST_P'; 'ST_F'; ...
                      'S2_TP'; 'S2_TF'; 'S2_PF'; ...
                      'S123'; 'Closure'; 'Max_S_minus_ST'};

        values_matrix = [S_all, ST_all, S2_all, S123_all, closure_all, violation_all];

        Trep.Mean = mean(values_matrix, 1)';
        Trep.Std = std(values_matrix, 0, 1)';
        Trep.Min = min(values_matrix, [], 1)';
        Trep.Max = max(values_matrix, [], 1)';

        disp(Trep);

        field_name = sprintf('t_%g', time);
        field_name = strrep(field_name, '.', '_');
        repeated_results.(field_name) = Trep;

        filename = sprintf('diagnostic_repeated_scrambling_t_%g.csv', time);
        writetable(Trep, filename);
        fprintf('Saved repeated scrambling table to %s\n', filename);

    end

    save('diagnostic_repeated_scrambling_results.mat', 'repeated_results');

end

%% ------------------------------------------------------------------------
% 5. OPTIONAL TIME-STEP SENSITIVITY TEST
% -------------------------------------------------------------------------

if run_sfe_tests && ~isempty(manual_time_steps_to_test)

    fprintf('\n============================================================\n');
    fprintf('5. TIME-STEP SENSITIVITY TEST\n');
    fprintf('============================================================\n');

    % Use final time point for time-step test by default
    time = test_time_points(end);

    % Use smaller N to avoid too much cost
    N_timestep = min(2^13, N_list(end));

    fprintf('Testing time-step sensitivity at t = %.0f min using N = %d\n', ...
        time, N_timestep);

    Tdt = table();

    for k_dt = 1:numel(manual_time_steps_to_test)

        time_step = manual_time_steps_to_test(k_dt);

        fprintf('\nRunning time_step = %.3g min...\n', time_step);

        cache = init_Extraction_Cache(time_step);
        model = @(YY) Simulate_Extraction_Cached(YY, time, time_step, cache);
        model_func = @(AA) applyToEachRowOptimized(model, AA);

        [A, B] = sobol_sampling_diagnostic(input_ranges, N_timestep, d, 8000 + k_dt);

        tic;
        [S, ST, S2, S123, YA, YB] = Sen_Saltelli( ...
            A, B, model_func, first_order_estimator, total_order_estimator);
        elapsed = toc;

        S2_pairs = get_second_order_pairs(S2);

        Tdt.time_step(k_dt,1) = time_step;
        Tdt.N(k_dt,1) = N_timestep;

        for i = 1:d
            Tdt.(sprintf('S_%s', NAMES{i}))(k_dt,1) = S(i);
            Tdt.(sprintf('ST_%s', NAMES{i}))(k_dt,1) = ST(i);
        end

        Tdt.S2_TP(k_dt,1) = S2(1,2);
        Tdt.S2_TF(k_dt,1) = S2(1,3);
        Tdt.S2_PF(k_dt,1) = S2(2,3);
        Tdt.S123(k_dt,1) = S123;
        Tdt.Closure(k_dt,1) = sum(S) + sum(S2_pairs) + S123;
        Tdt.Max_S_minus_ST(k_dt,1) = max(S - ST);
        Tdt.OutputMean(k_dt,1) = mean([YA; YB]);
        Tdt.OutputStd(k_dt,1) = std([YA; YB]);
        Tdt.Elapsed_s(k_dt,1) = elapsed;

    end

    disp(Tdt);

    writetable(Tdt, sprintf('diagnostic_timestep_sensitivity_t_%g.csv', time));
    save('diagnostic_timestep_sensitivity_results.mat', 'Tdt');

end

fprintf('\n============================================================\n');
fprintf('DIAGNOSTICS COMPLETE\n');
fprintf('============================================================\n');

fprintf('\nRecommended checks before full rerun:\n');
fprintf('  1. Benchmark additive model: S_i should be close to ST_i.\n');
fprintf('  2. Benchmark additive model: S2 and S123 should be near zero.\n');
fprintf('  3. In SFE model, large positive max(S-ST) indicates a problem.\n');
fprintf('  4. Important indices should stabilize as N increases.\n');
fprintf('  5. Repeated scrambled sequences should give similar rankings.\n');
fprintf('  6. Time-step changes should not strongly change ST indices.\n');

%% ========================================================================
% LOCAL HELPER FUNCTIONS
% ========================================================================

function [A, B] = sobol_sampling_diagnostic(input_ranges, N, d, seed)
% Generate scrambled Sobol matrices A and B.
%
% This version uses a variable seed through rng before scrambling.
% Different seeds give independent randomized QMC replicates.

    rng(seed, 'twister');

    sob = sobolset(2*d, 'Skip', 1000, 'Leap', 100);
    sob = scramble(sob, 'MatousekAffineOwen');

    samples = net(sob, N);

    A01 = samples(:, 1:d);
    B01 = samples(:, d+1:2*d);

    A = scale_to_ranges_diagnostic(A01, input_ranges);
    B = scale_to_ranges_diagnostic(B01, input_ranges);
end

function X = scale_to_ranges_diagnostic(X01, ranges)
% Scale [0,1] samples to physical input ranges.

    mins = ranges(:,1)';
    maxs = ranges(:,2)';
    X = X01 .* (maxs - mins) + mins;
end

function time_step = select_time_steps_diagnostic(t)
% Same time-step logic as main script, but for one scalar time.

    if t <= 30
        time_step = 0.5;
    elseif t <= 150
        time_step = 1;
    elseif t <= 300
        time_step = 2;
    elseif t <= 600
        time_step = 5;
    elseif t <= 1000
        time_step = 10;
    else
        time_step = 15;
    end
end

function S2_pairs = get_second_order_pairs(S2)
% Extract unique second-order terms for d = 3.

    S2_pairs = [S2(1,2), S2(1,3), S2(2,3)];
end