function [S, ST, estimators] = Sen_Sobol_Estimators(A, B, f, varargin)
% SEN_SOBOL_ESTIMATORS  Compute Sobol indices using multiple estimators
%
%   [S, ST, estimators] = Sen_Sobol_Estimators(A, B, f, ...)
%
%   Computes first-order (S) and total-effect (ST) Sobol indices using
%   three different estimator methods for comparison.
%
%   Inputs:
%       A, B : N-by-k sample matrices (same size). Quasi-random recommended.
%       f    : function handle; must return N-by-1 vector for f(X,...)
%       ...  : extra args passed to f via varargin
%
%   Outputs:
%       S    : 1-by-k vector of first-order indices (Janon-Monod estimator)
%       ST   : 1-by-k vector of total-effect indices (Jansen estimator)
%       estimators : struct with all estimator results for comparison
%           .S_Saltelli   : First-order (Saltelli 2010)
%           .S_Janon      : First-order (Janon-Monod 2014) - recommended
%           .S_Sobol      : First-order (Sobol 1993, original)
%           .ST_Jansen    : Total-effect (Jansen 1999) - recommended
%           .ST_Sobol     : Total-effect (Sobol 2001)
%           .ST_Saltelli  : Total-effect (Saltelli 2008)
%           .yA, .yB      : Model outputs for samples A and B
%           .VarY         : Total output variance
%           .f0           : Mean output
%
%   References:
%       [1] Saltelli et al. (2010) "Variance based sensitivity analysis"
%       [2] Janon et al. (2014) "Asymptotic normality and efficiency"
%       [3] Jansen (1999) "Analysis of variance designs for model output"
%       [4] Sobol (1993) "Sensitivity estimates for nonlinear models"
%
%   Note: Janon-Monod (S) and Jansen (ST) are numerically most stable.

    [N, k] = size(A);
    assert(isequal(size(B), [N, k]), 'A and B must be N-by-k.');

    %% Evaluate base matrices
    fprintf('Evaluating base matrices A and B (%d samples each)...\n', N);
    tic
    yA = f(A, varargin{:}); yA = yA(:);
    yB = f(B, varargin{:}); yB = yB(:);
    fprintf('Base evaluations completed in %.2f seconds\n', toc);

    %% Compute variance statistics
    yAll = [yA; yB];
    f0 = mean(yAll);
    VarY = var(yAll);  % Using MATLAB var (N-1 denominator)

    % Alternative variance (population, N denominator)
    VarY_pop = mean(yAll.^2) - f0^2;

    if VarY <= 0
        error('Estimated output variance <= 0. Check model outputs or sample size.');
    end

    %% Initialize storage for all estimators
    S_Saltelli = zeros(1, k);
    S_Janon = zeros(1, k);
    S_Sobol = zeros(1, k);

    ST_Jansen = zeros(1, k);
    ST_Sobol = zeros(1, k);
    ST_Saltelli = zeros(1, k);

    %% Compute indices for each parameter
    fprintf('Computing Sobol indices for %d parameters...\n', k);

    for i = 1:k
        tic

        % Hybrid matrix A_Bi: A with column i replaced by B's column i
        ABi = A; ABi(:, i) = B(:, i);
        yABi = f(ABi, varargin{:}); yABi = yABi(:);

        % ================================================================
        % FIRST-ORDER ESTIMATORS (S_i)
        % ================================================================

        % 1. Saltelli (2010) estimator
        %    S_i = (1/N * sum(yB .* yABi) - f0^2) / Var(Y)
        S_Saltelli(i) = (mean(yB .* yABi) - f0^2) / VarY_pop;

        % 2. Janon-Monod (2014) estimator - RECOMMENDED
        %    S_i = (1/N * sum(yB .* (yABi - yA))) / Var(Y)
        %    More stable: avoids subtracting f0^2 (cancellation errors)
        S_Janon(i) = mean(yB .* (yABi - yA)) / VarY;

        % 3. Sobol (1993) original estimator
        %    S_i = (1/N * sum(yA .* (yABi - yB))) / Var(Y)
        %    Note: Uses different pairing than Janon
        S_Sobol(i) = mean(yA .* (yABi - yB)) / VarY;

        % ================================================================
        % TOTAL-EFFECT ESTIMATORS (ST_i)
        % ================================================================

        % 1. Jansen (1999) estimator - RECOMMENDED
        %    ST_i = (1/2N * sum((yA - yABi).^2)) / Var(Y)
        %    Most stable: uses squared differences
        ST_Jansen(i) = mean((yA - yABi).^2) / (2 * VarY);

        % 2. Sobol (2001) estimator
        %    ST_i = 1 - (1/N * sum(yB .* yABi) - f0^2) / Var(Y)
        %    Note: ST = 1 - S_~i (complementary first-order)
        ST_Sobol(i) = 1 - (mean(yB .* yABi) - f0^2) / VarY_pop;

        % 3. Saltelli (2008) estimator
        %    ST_i = (1/N * sum(yA .* (yA - yABi))) / Var(Y)
        ST_Saltelli(i) = mean(yA .* (yA - yABi)) / VarY;

        fprintf('  Parameter %d/%d completed (%.2f sec)\n', i, k, toc);
    end

    %% Return recommended estimators as primary output
    S = S_Janon;      % Janon-Monod is most stable for first-order
    ST = ST_Jansen;   % Jansen is most stable for total-effect

    %% Store all estimators for comparison
    estimators.S_Saltelli = S_Saltelli;
    estimators.S_Janon = S_Janon;
    estimators.S_Sobol = S_Sobol;

    estimators.ST_Jansen = ST_Jansen;
    estimators.ST_Sobol = ST_Sobol;
    estimators.ST_Saltelli = ST_Saltelli;

    estimators.yA = yA;
    estimators.yB = yB;
    estimators.VarY = VarY;
    estimators.f0 = f0;

    %% Print comparison summary
    fprintf('\n--- Estimator Comparison ---\n');
    fprintf('First-order (S):\n');
    fprintf('  Saltelli:   [%s]\n', sprintf('%.4f ', S_Saltelli));
    fprintf('  Janon:      [%s] (recommended)\n', sprintf('%.4f ', S_Janon));
    fprintf('  Sobol:      [%s]\n', sprintf('%.4f ', S_Sobol));
    fprintf('Total-effect (ST):\n');
    fprintf('  Jansen:     [%s] (recommended)\n', sprintf('%.4f ', ST_Jansen));
    fprintf('  Sobol:      [%s]\n', sprintf('%.4f ', ST_Sobol));
    fprintf('  Saltelli:   [%s]\n', sprintf('%.4f ', ST_Saltelli));
    fprintf('Interaction (ST - S) using Jansen/Janon:\n');
    fprintf('  [%s]\n', sprintf('%.4f ', ST_Jansen - S_Janon));

end
