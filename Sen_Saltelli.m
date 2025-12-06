function [S, ST, yA, yB] = Sen_Saltelli(A, B, f, varargin)
% SEN_SALTELLI  First-order (S) and total-effect (ST) Sobol' indices.
% Inputs:
%   A, B : N-by-k sample matrices (same size). Quasi-random recommended.
%   f    : function handle; must return N-by-1 (or 1-by-N) vector for f(X,...)
%   ...  : extra args passed to f via varargin
% Outputs:
%   S    : 1-by-k vector of first-order indices
%   ST   : 1-by-k vector of total-effect indices
%
% Estimators:
%   - S  : Saltelli single-loop with hybrids C_i = B; C_i(:,i)=A(:,i), paired with yA
%   - ST : Jansen (1999) estimator (numerically stable)

    [N, k] = size(A);
    assert(isequal(size(B), [N, k]), 'A and B must be N-by-k.');

    tic
    yA = f(A, varargin{:}); yA = yA(:);
    yB = f(B, varargin{:}); yB = yB(:);
    toc
    
    yAll = [yA; yB];
    f0   = mean(yAll);
    VarY = mean(yAll.^2) - f0^2;
    if VarY <= 0
        error('Estimated output variance <= 0. Check model outputs or sample size.');
    end

    S  = zeros(1, k);
    ST = zeros(1, k);

    for i = 1:k
        Ci = B; Ci(:, i) = A(:, i);             % hybrid B_Ai
        yCi = f(Ci, varargin{:}); yCi = yCi(:);

        % First-order (Saltelli single-loop)
        S(i)  = (mean(yA .* yCi) - f0^2) / VarY;

        % Total-effect (Jansen)
        ST(i) = mean((yB - yCi).^2) / (2 * VarY);

        % Optional progress:
        %fprintf('%d/%d: S=%.6g, ST=%.6g\n', i, k, S(i), ST(i));
    end
end
