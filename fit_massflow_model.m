function results = fit_massflow_model(Y, F, mu, D, A)
% FIT_MASSFLOW_MODEL
%   Fits:  Y = Y0 + F*(Y1 + Y2*D/(A*mu))
%   Unknowns: [Y0; Y1; Y2]
%   Known per-row: F(i), mu(i); constants: D, A
%
%   results fields:
%     .beta        [Y0; Y1; Y2] estimates
%     .cov         variance-covariance matrix of beta
%     .se          standard errors of beta
%     .tstat, .p   t-stats and p-values (dof = n-3)
%     .ci95        95% confidence intervals for beta
%     .goodness    struct with R2, adjR2, RMSE, SSE, SSR, SST, F, pF, dof
%     .design      design matrix used (X)
%
% Example:
%   results = fit_massflow_model(Y, F, mu, D, A);

    %--- basic checks & clean ---
    Y  = Y(:); F = F(:); mu = mu(:);
    ok = isfinite(Y) & isfinite(F) & isfinite(mu);
    Y = Y(ok); F = F(ok); mu = mu(ok);

    if ~isscalar(D) || ~isscalar(A)
        error('D and A must be scalars.');
    end

    n = numel(Y);
    if n < 3
        error('Need at least 3 observations.');
    end

    %--- design matrix for linear least squares ---
    % Y = [1, F, F.*(D/(A*mu))] * [Y0; Y1; Y2]
    k = D / A;
    X = [ones(n,1), F, F .* (k ./ mu)];
    p = size(X,2);

    %--- OLS fit ---
    beta = X \ Y;                         % estimates
    yhat = X * beta;
    r    = Y - yhat;

    %--- sums of squares & GOF ---
    SSE  = sum(r.^2);
    ybar = mean(Y);
    SST  = sum((Y - ybar).^2);
    SSR  = SST - SSE;

    dof  = n - p;
    s2   = SSE / dof;                     % residual variance
    RMSE = sqrt(s2);
    R2   = 1 - SSE / SST;
    adjR2 = 1 - (1 - R2) * (n - 1) / (n - p);

    % Overall F-test for regression (with intercept)
    dfr = p - 1;
    Fstat = (SSR / dfr) / (SSE / dof);
    pF    = 1 - fcdf(Fstat, dfr, dof);

    %--- variance-covariance of parameters & SEs ---
    XtX = X' * X;
    CovBeta = s2 * (XtX \ eye(p));        % numerically stable inverse
    se = sqrt(diag(CovBeta));

    %--- inference per parameter ---
    tstat = beta ./ se;
    pval  = 2 * (1 - tcdf(abs(tstat), dof));
    tcrit = tinv(0.975, dof);
    ci95  = [beta - tcrit*se, beta + tcrit*se];

    %--- pack results ---
    names = {'Y0','Y1','Y2'}';
    results.beta     = beta;
    results.se       = se;
    results.cov      = CovBeta;
    results.tstat    = tstat;
    results.p        = pval;
    results.ci95     = ci95;
    results.design   = X;
    results.goodness = struct('n',n,'p',p,'dof',dof,'SSE',SSE,'SSR',SSR,'SST',SST, ...
                              'RMSE',RMSE,'R2',R2,'AdjR2',adjR2,'F',Fstat,'pF',pF);
    results.paramTable = table(beta,se,tstat,pval,ci95(:,1),ci95(:,2), ...
                               'VariableNames',{'Estimate','SE','tStat','pValue','CI95_L','CI95_U'}, ...
                               'RowNames',names);

    %--- display summary ---
    fprintf('\n=== Fit summary: Y = Y0 + F*(Y1 + Y2*D/(A*mu)) ===\n');
    disp(results.paramTable);
    fprintf('RMSE = %.6g,   R^2 = %.6f,   Adj R^2 = %.6f\n', RMSE, R2, adjR2);
    fprintf('SSE = %.6g, SSR = %.6g, F(%d,%d) = %.3f, p = %.3g\n', ...
            SSE, SSR, dfr, dof, Fstat, pF);
    fprintf('\nVariance–Covariance matrix of [Y0 Y1 Y2]^T:\n');
    disp(array2table(CovBeta,'VariableNames',names,'RowNames',names));
end
