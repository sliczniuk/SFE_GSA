function [S, ST, S2, S3, yA, yB, CI, evals] = Sen_Saltelli_TimeSeries(A, B, f, method_s, method_st, varargin)
% SEN_SALTELLI_TIMESERIES Sobol indices for vector-valued time outputs.
%   f must return an N-by-nTimes matrix when called with an N-by-k sample
%   matrix. Each column is treated as one time point.

[N, k] = size(A);
assert(isequal(size(B), [N, k]), 'A and B must be N-by-k.');

n_bootstrap = 0;
confidence_level = 0.95;

idx = 1;
while idx <= numel(varargin)
    if ischar(varargin{idx}) && strcmpi(varargin{idx}, 'Bootstrap')
        if idx + 1 > numel(varargin) || ~isnumeric(varargin{idx + 1})
            error('Bootstrap option must be followed by a numeric resample count.');
        end
        n_bootstrap = varargin{idx + 1};
        varargin(idx:idx + 1) = [];
    else
        idx = idx + 1;
    end
end

n_hybrid_i = k;
n_hybrid_ij = k * (k - 1) / 2;
evals.total_model_rows = (2 + n_hybrid_i + n_hybrid_ij) * N;
evals.base_rows = 2 * N;
evals.first_total_rows = n_hybrid_i * N;
evals.second_order_rows = n_hybrid_ij * N;

fprintf('\n=== TIME-SERIES SOBOL SENSITIVITY ANALYSIS ===\n');
fprintf('Sample size N = %d, Parameters k = %d\n', N, k);
fprintf('Model row evaluations needed: %d\n', evals.total_model_rows);

fprintf('\n[1/3] Evaluating base matrices A and B...\n');
tic_base = tic;
yA = f(A, varargin{:});
yB = f(B, varargin{:});
yA = ensure_2d_output(yA, N, 'A');
yB = ensure_2d_output(yB, N, 'B');
n_times = size(yA, 2);
if size(yB, 2) ~= n_times
    error('Model outputs for A and B must have the same number of time columns.');
end
base_time = toc(tic_base);
fprintf('  Completed in %.2f seconds (%.0f rows/sec)\n', base_time, 2 * N / base_time);

if any(~isfinite(yA(:)))
    warning('Matrix A produced %d non-finite outputs (NaN/Inf).', sum(~isfinite(yA(:))));
end
if any(~isfinite(yB(:)))
    warning('Matrix B produced %d non-finite outputs (NaN/Inf).', sum(~isfinite(yB(:))));
end

yABi_all = zeros(N, n_times, k);

fprintf('\n[2/3] Evaluating A_Bi matrices for first/total indices...\n');
tic_first = tic;
for i = 1:k
    fprintf('  [%d/%d] Processing parameter %d...', i, k, i);
    ABi = A;
    ABi(:, i) = B(:, i);
    yABi = ensure_2d_output(f(ABi, varargin{:}), N, sprintf('A_B%d', i));
    if size(yABi, 2) ~= n_times
        error('Model output for A_B%d has %d columns; expected %d.', i, size(yABi, 2), n_times);
    end
    yABi_all(:, :, i) = yABi;
    fprintf(' done\n');
end
fprintf('  First/total model outputs completed in %.2f seconds\n', toc(tic_first));

yABij_all = NaN(N, n_times, k, k);

fprintf('\n[3/3] Evaluating A_Bij matrices for second-order indices...\n');
tic_second = tic;
pair_count = 0;
for i = 1:(k - 1)
    for j = (i + 1):k
        pair_count = pair_count + 1;
        fprintf('  [%d/%d] Processing pair (%d,%d)...', pair_count, n_hybrid_ij, i, j);
        ABij = A;
        ABij(:, i) = B(:, i);
        ABij(:, j) = B(:, j);
        yABij = ensure_2d_output(f(ABij, varargin{:}), N, sprintf('A_B%d%d', i, j));
        if size(yABij, 2) ~= n_times
            error('Model output for A_B%d%d has %d columns; expected %d.', i, j, size(yABij, 2), n_times);
        end
        yABij_all(:, :, i, j) = yABij;
        fprintf(' done\n');
    end
end
fprintf('  Second-order model outputs completed in %.2f seconds\n', toc(tic_second));

[S, ST, S2, S3] = compute_indices(yA, yB, yABi_all, yABij_all, method_s, method_st);

CI = empty_ci(n_times, k, confidence_level, n_bootstrap);
if n_bootstrap > 0
    fprintf('\n=== BOOTSTRAP CONFIDENCE INTERVALS ===\n');
    fprintf('Resampling %d times for %.0f%% CI...\n', n_bootstrap, confidence_level * 100);

    S_boot = NaN(n_bootstrap, n_times, k);
    ST_boot = NaN(n_bootstrap, n_times, k);
    S2_boot = NaN(n_bootstrap, n_times, k, k);
    S3_boot = NaN(n_bootstrap, n_times);

    tic_boot = tic;
    for b = 1:n_bootstrap
        boot_idx = randi(N, N, 1);
        try
            [S_b, ST_b, S2_b, S3_b] = compute_indices( ...
                yA(boot_idx, :), ...
                yB(boot_idx, :), ...
                yABi_all(boot_idx, :, :), ...
                yABij_all(boot_idx, :, :, :), ...
                method_s, method_st);
        catch ME
            if contains(ME.message, 'Estimated output variance')
                continue;
            end
            rethrow(ME);
        end

        S_boot(b, :, :) = reshape(S_b, [1, n_times, k]);
        ST_boot(b, :, :) = reshape(ST_b, [1, n_times, k]);
        S2_boot(b, :, :, :) = reshape(S2_b, [1, n_times, k, k]);
        S3_boot(b, :) = reshape(S3_b, [1, n_times]);

        if mod(b, max(1, floor(n_bootstrap / 10))) == 0
            elapsed = toc(tic_boot);
            remaining = (n_bootstrap - b) * elapsed / b;
            fprintf('  Progress: %d/%d (%.1f%%) - ETA: %.1f sec\n', ...
                b, n_bootstrap, 100 * b / n_bootstrap, remaining);
        end
    end

    alpha = 1 - confidence_level;
    pct = [100 * alpha / 2, 100 * (1 - alpha / 2)];

    for t = 1:n_times
        for i = 1:k
            CI.S_CI(t, i, :) = percentile_finite(S_boot(:, t, i), pct);
            CI.ST_CI(t, i, :) = percentile_finite(ST_boot(:, t, i), pct);
        end
        for i = 1:(k - 1)
            for j = (i + 1):k
                CI.S2_CI(t, i, j, :) = percentile_finite(S2_boot(:, t, i, j), pct);
            end
        end
        if k == 3
            CI.S3_CI(t, :) = percentile_finite(S3_boot(:, t), pct);
        end
    end

    fprintf('  Bootstrap completed in %.2f seconds\n', toc(tic_boot));
end

end

function y = ensure_2d_output(y, expected_rows, label)
y = double(y);
if isvector(y)
    y = y(:);
end
if size(y, 1) ~= expected_rows
    error('Model output for %s has %d rows; expected %d.', label, size(y, 1), expected_rows);
end
end

function [S, ST, S2, S3] = compute_indices(yA, yB, yABi_all, yABij_all, method_s, method_st)
[N, n_times] = size(yA); %#ok<ASGLU>
k = size(yABi_all, 3);

yAll = [yA; yB];
f0 = mean(yAll, 1);
VarY = mean((yAll - f0) .^ 2, 1);
if any(VarY <= 1e-12)
    bad_time = find(VarY <= 1e-12, 1, 'first');
    error('Estimated output variance <= 1e-12 at time column %d (%.6e).', bad_time, VarY(bad_time));
end

S = zeros(n_times, k);
ST = zeros(n_times, k);
S2 = NaN(n_times, k, k);
S3 = NaN(n_times, 1);

for i = 1:k
    yABi = yABi_all(:, :, i);

    if strcmpi(method_s, 'Saltelli')
        S(:, i) = ((mean(yB .* yABi, 1) - f0 .^ 2) ./ VarY)';
    elseif strcmpi(method_s, 'Janon-Monod')
        S(:, i) = (mean(yB .* (yABi - yA), 1) ./ VarY)';
    elseif strcmpi(method_s, 'Sobol')
        S(:, i) = (mean(yA .* (yABi - yB), 1) ./ VarY)';
    else
        error('Incorrect method_s specified: "%s". Use Saltelli, Janon-Monod, or Sobol.', method_s);
    end

    if strcmpi(method_st, 'Jansen')
        ST(:, i) = (mean((yA - yABi) .^ 2, 1) ./ (2 * VarY))';
    elseif strcmpi(method_st, 'Sobol')
        ST(:, i) = (1 - (mean(yB .* yABi, 1) - f0 .^ 2) ./ VarY)';
    elseif strcmpi(method_st, 'Saltelli')
        ST(:, i) = (mean(yA .* (yA - yABi), 1) ./ VarY)';
    else
        error('Incorrect method_st specified: "%s". Use Jansen, Sobol, or Saltelli.', method_st);
    end
end

for i = 1:(k - 1)
    for j = (i + 1):k
        yABij = yABij_all(:, :, i, j);
        if strcmpi(method_s, 'Janon-Monod')
            S2(:, i, j) = (mean(yB .* (yABij - yABi_all(:, :, i) - yABi_all(:, :, j) + yA), 1) ./ VarY)';
        else
            closed_second = (mean(yB .* yABij, 1) - f0 .^ 2) ./ VarY;
            S2(:, i, j) = closed_second' - S(:, i) - S(:, j);
        end
    end
end

if k == 3
    yAB1 = yABi_all(:, :, 1);
    yAB2 = yABi_all(:, :, 2);
    yAB3 = yABi_all(:, :, 3);
    yAB12 = yABij_all(:, :, 1, 2);
    yAB13 = yABij_all(:, :, 1, 3);
    yAB23 = yABij_all(:, :, 2, 3);

    if strcmpi(method_s, 'Janon-Monod')
        % Direct third-order inclusion-exclusion estimate. A_B123 equals B,
        % so yB is reused and no extra model evaluations are needed.
        S3 = (mean(yB .* (yB - yAB12 - yAB13 - yAB23 + ...
            yAB1 + yAB2 + yAB3 - yA), 1) ./ VarY)';
    else
        closed_all = (mean(yB .* yB, 1) - f0 .^ 2) ./ VarY;
        for t = 1:n_times
            S2_t = squeeze(S2(t, :, :));
            S3(t) = closed_all(t) - sum(S(t, :)) - sum(S2_t(~isnan(S2_t)));
        end
    end
end
end

function CI = empty_ci(n_times, k, confidence_level, n_bootstrap)
CI.S_CI = [];
CI.ST_CI = [];
CI.S2_CI = [];
CI.S3_CI = [];
CI.confidence_level = confidence_level;
CI.n_bootstrap = n_bootstrap;

if n_bootstrap > 0
    CI.S_CI = NaN(n_times, k, 2);
    CI.ST_CI = NaN(n_times, k, 2);
    CI.S2_CI = NaN(n_times, k, k, 2);
    CI.S3_CI = NaN(n_times, 2);
end
end

function ci = percentile_finite(values, pct)
values = values(isfinite(values));
if isempty(values)
    ci = [NaN, NaN];
else
    ci = prctile(values, pct);
end
end
