% Sobol Sensitivity Analysis in MATLAB - OPTIMIZED VERSION
% Major performance improvements:
%   1. Parameters.csv loaded ONCE (not 50,000 times)
%   2. CasADi integrator built ONCE per time point (not per sample)
%   3. Sobol quasi-random sequences for better convergence
%   4. Preallocated numeric arrays instead of cells

close all; clear; clc;

%% Initialize parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

%% USER INPUTS
alpha = 0.5;
FONT = 10;

% Number of Monte Carlo samples (can use fewer with quasi-random)
N = 1e4;

NAMES = {'T[$^\circ C$]', 'P[bar]', 'F[kg/s]'};

input_ranges = [30   40;
                100  200;
                3.33 6.67];

% Number of input parameters
d = size(input_ranges, 1);

%% Generate Sobol quasi-random samples (better convergence than pseudo-random)
[A, B] = sobol_sampling(input_ranges, N);

diary('myTextLog_optimized.txt');

%% Time points to analyze
time_points = [5, 15, 30, 60, 90, 120, 150, 240, 300, 450, 600, 750, 900, 1200, 1500, 2000];

for idx = 1:length(time_points)
    time = time_points(idx);

    fprintf('\n========================================\n');
    fprintf('Processing TIME = %.0f [min]\n', time);
    fprintf('========================================\n');

    %% Initialize cache for this time point (loads params & builds integrator ONCE)
    tic
    cache = init_Extraction_Cache(1);  % time_step = 1 minute
    cache_time = toc;
    fprintf('Cache initialization: %.2f seconds\n', cache_time);

    %% Create model function using cached data
    model = @(YY) Simulate_Extraction_Cached(YY, time, 1, cache);
    model_func = @(A) applyToEachRowOptimized(model, A);

    %% Run Sobol analysis
    tic
    [first_order, total_order, YA, YB] = Sen_Saltelli(A, B, model_func);
    Y = [YA; YB];
    VY = var(Y);
    analysis_time = toc;
    fprintf('Sobol analysis: %.2f seconds\n', analysis_time);

    %% Display results
    fprintf('\nTIME = %.0f [min]\n', time);
    fprintf('\nSum of First-order indices: %.4f\n', sum(first_order));
    fprintf('Sum of Total-order indices: %.4f\n\n', sum(total_order));

    % Interaction terms (should be non-negative)
    interaction = total_order - first_order;
    fprintf('Interaction terms: [%.4f, %.4f, %.4f]\n', interaction(1), interaction(2), interaction(3));

    fprintf('\nSobol Sensitivity Indices:\n');
    for i = 1:d
        fprintf('%s: First-order = %.4f, Total-order = %.4f, Interaction = %.4f\n', ...
            NAMES{i}, first_order(i), total_order(i), interaction(i));
    end

    %% Scatter plots
    for i = 1:d
        x = A(:,i);
        y = YA;

        coeffs = polyfit(A(:,i), YA, 1);
        xfit = [min(x) max(x)];
        yfit = polyval(coeffs, xfit);

        yCalc1 = polyval(coeffs, x);
        Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);

        if i == 3
            xfit = xfit * 1e-5;
            x = x * 1e-5;
        end

        figure(1);
        hold on
        scatter(x, y, 10, y, 'filled', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
        colorbar;
        plot(xfit, yfit, 'k-', 'LineWidth', 2);
        hold off

        title(sprintf('Scatter plot after %.0f [min]\n $y = %.6f \\cdot x + %.6f, R^2 = %.2f$', time, coeffs, Rsq1))
        xlabel(sprintf('%s', NAMES{i}));
        ylabel('Yield [g]');

        grid off; axis square;
        set(gca, 'FontSize', FONT)
        exportgraphics(figure(1), ['GSA_Scatter_' + string(NAMES{i}(1)) + '_' + string(time) + '.png'], 'Resolution', 300);
        close all;
    end

    %% Output distribution
    all_outputs = [YA; YB];

    figure(1);
    histogram(all_outputs, 'Normalization', 'pdf', 'FaceColor', [0.2, 0.6, 0.8]);
    xlabel('$Yield [g]$');
    ylabel('Probability Density');
    title(sprintf('Probability density plot after %.0f [min]', time))
    grid off;

    output_mean = mean(all_outputs);
    output_std = std(all_outputs);

    fprintf('\nOutput Statistics:\n');
    fprintf('Mean of output: %.4f\n', output_mean);
    fprintf('Standard deviation of output: %.4f\n', output_std);
    set(gca, 'FontSize', FONT)
    exportgraphics(figure(1), ['GSA_Distribution_' + string(time) + '.png'], 'Resolution', 300);
    close all;

    fprintf('\nTotal time for TIME=%.0f: %.2f seconds\n', time, cache_time + analysis_time);
end

diary('off');
fprintf('\n\nAnalysis complete!\n');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

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
