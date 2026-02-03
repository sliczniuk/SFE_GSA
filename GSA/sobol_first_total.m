function [first_order, total_order] = sobol_first_total(model_func, A, B)
% Computes first-order and total-order Sobol indices using Saltelli's method (no parfor)
% 
% Inputs:
%   model_func - function handle: y = model_func(x), where x is a 1xd vector
%   A, B       - Nxd input sample matrices
%
% Outputs:
%   first_order    - dx1 vector of first-order Sobol indices
%   total_order    - dx1 vector of total-order Sobol indices

N = size(A, 1);
d = size(A, 2);

% Preallocate
YA  = zeros(N, 1);
YB  = zeros(N, 1);
YAB = zeros(N, d);
first_order = zeros(d,1);
total_order = zeros(d,1);

% Evaluate model at A and B
parfor i = 1:N
    YA(i) = model_func(A(i, :));
    YB(i) = model_func(B(i, :));
end

Y = [YA; YB];
f0 = mean(Y);
VY = var(Y);

%% STEP 2: Compute First-order and Total-order Sobol indices, store YABi

for i = 1:d
    % Hybrid matrix A_Bi: A with column i replaced by B's column i
    ABi = A;
    ABi(:, i) = B(:, i);
    parfor k = 1:N
        YAB(k,i) = model_func(ABi(k,:));
    end

    % First-order (Saltelli 2010): pair yB with y(A_Bi)
    first_order(i) = (mean(YB .* YAB(:,i)) - f0^2) / VY;

    % Total-effect (Jansen 1999): pair yA with y(A_Bi)
    total_order(i) = mean((YA - YAB(:,i)).^2) / (2 * VY);
end

end
