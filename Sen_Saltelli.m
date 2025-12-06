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
% Estimators (Saltelli 2010):
%   - S  : First-order using hybrid A_Bi (A with col i from B), paired with yB
%   - ST : Total-effect (Jansen 1999), using (yA - y(A_Bi))^2

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
        % Hybrid matrix A_Bi: A with column i replaced by B's column i
        ABi = A; ABi(:, i) = B(:, i);
        yABi = f(ABi, varargin{:}); yABi = yABi(:);

        % First-order (Saltelli 2010, Eq. in "Variance based sensitivity analysis")
        % S_i = V[E[Y|X_i]] / V[Y]
        S(i)  = (mean(yB .* yABi) - f0^2) / VarY;

        % Total-effect (Jansen 1999 estimator)
        % ST_i = E[(Y(A) - Y(A_Bi))^2] / (2*V[Y])
        ST(i) = mean((yA - yABi).^2) / (2 * VarY);

        % Optional progress:
        %fprintf('%d/%d: S=%.6g, ST=%.6g\n', i, k, S(i), ST(i));
    end
end
