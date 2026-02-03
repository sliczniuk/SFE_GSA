% Sobol Sensitivity Analysis in MATLAB
% General-purpose script applicable to any deterministic model

close all; clear; clc;

if isempty(gcp('nocreate'))
    parpool;
end

%% USER INPUTS

alpha      = 0.5;
FONT       = 10;

% Number of Monte Carlo samples (recommended: thousands for accuracy)
N = 1e4;

NAMES     = {'T[$^\circ C$]', 'P[bar]', 'F[kg/s]'};

input_ranges = [30   40;
                100  200;
                3.33 6.67];
    
[A, B] = saltelli_sampling(input_ranges, N);

% Number of input parameters (dimensions)
d = size(input_ranges,1);

diary('myTextLog.txt');

delta_time = [0.01, 0.1,  1, 1 ,  1,   1,   1,  5,    5,  10,  10,  10,  15,  30,   25,   50];

for time = [  5   , 15 , 30, 60, 90, 120, 150, 240, 300, 450, 600, 750, 900, 1200, 1500, 2000]

    model      = @(YY) Simulate_Extraction(YY, time, 1);
    model_func = @(A) applyToEachRow(model, A);
    %model_func = @(A) applyToEachRowGPU_arrayfun(model, A);
    %model_func = @(A) applyToEachRowGPU_arrayfun(time, A);
    
    tic
    [first_order, total_order, YA, YB] = Sen_Saltelli(A, B, model_func);
    Y = [YA; YB];
    VY = var(Y);
    toc
    
    %%
    fprintf('\nTIME = %.0f [min]\n', time);
    fprintf('\nSum of First-order indices: %.4f\n', sum(first_order));
    fprintf('Sum of Total-order indices: %.4f\n\n', sum(total_order));
    
    %% STEP 3: Display Results
    
    fprintf('\nSobol Sensitivity Indices:\n');
    for i = 1:d
        fprintf('%s: First-order = %.4f, Total-order = %.4f\n', ...
            NAMES{i}, first_order(i), total_order(i));
    end
    
    for i = 1:d
        
        x = A(:,i);
        y = YA;
        
        coeffs = polyfit(A(:,i), YA, 1);
        xfit = [min(x) max(x)];
        yfit = polyval(coeffs, xfit);
    
        yCalc1 = polyval(coeffs, x);
        Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2);
        
        if i==3
            xfit = xfit * 1e-5;
            x    = x    * 1e-5;
        end

        figure(1);
        hold on
        scatter(x, y, 10, y, 'filled', 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha); colorbar;
        plot(xfit, yfit, 'k-', LineWidth=2);
        hold off

        title(sprintf('Scatter plot after %.0f [min]\n $y = %.6f \\cdot x + %.6f, R^2 = %.2f$', time, coeffs, Rsq1))
        xlabel(sprintf('%s', NAMES{i}));
        ylabel('Yield [g]');
        
        grid off; axis square;
        set(gca,'FontSize',FONT)
        exportgraphics(figure(1), ['GSA_Scatter_'+string(NAMES{i}(1))+'_'+string(time)+'.png'], "Resolution",300); close all;
    end
    
    %% STEP 7: Compute and Plot Output Distribution
    
    % Compute output samples (already computed YA and YB can be combined)
    all_outputs = [YA; YB];
    
    % Plot histogram of output distribution
    figure(1);
    histogram(all_outputs, 'Normalization', 'pdf', 'FaceColor', [0.2, 0.6, 0.8]);
    xlabel('$Yield [g]$');
    ylabel('Probability Density');
    title(sprintf('Probability density plot after %.0f [min]', time))
    grid off;
    
    % Compute statistics (mean and standard deviation)
    output_mean = mean(all_outputs);
    output_std = std(all_outputs);
    
    % Display statistics
    fprintf('\nOutput Statistics:\n');
    fprintf('Mean of output: %.4f\n', output_mean);
    fprintf('Standard deviation of output: %.4f\n', output_std);
    set(gca,'FontSize',FONT)
    exportgraphics(figure(1), ['GSA_Distribution_'+string(time)+'.png'], "Resolution",300); close all;
    
end

diary('off');

%%

function [A, B] = saltelli_sampling(input_ranges, N)
    % SALTELLI_SAMPLING Generate matrices A and B for Sobol-Saltelli GSA
    %
    %   INPUTS:
    %       input_ranges : k-by-2 matrix, ranges for each variable [min, max]
    %       N            : number of Monte Carlo samples (rows)
    %
    %   OUTPUTS:
    %       A, B         : N-by-k matrices of uniform samples
    %
    % Example:
    %   ranges = [0 1;    % x1 in [0,1]
    %             -2 3;   % x2 in [-2,3]
    %             10 20]; % x3 in [10,20]
    %   [A, B] = saltelli_sampling(ranges, 1000);
    
        % Number of input variables
        k = size(input_ranges, 1);
    
        % Random samples in [0,1]
        A01 = rand(N, k);
        B01 = rand(N, k);
    
        % Scale to desired ranges
        A = scale_to_ranges(A01, input_ranges);
        B = scale_to_ranges(B01, input_ranges);
    end
    
    
    function X = scale_to_ranges(X01, ranges)
    % Helper: scale a [0,1] sample to user-defined ranges
        mins = ranges(:,1)';
        maxs = ranges(:,2)';
        X = X01 .* (maxs - mins) + mins;
    end