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
VY = var(Y);

%% STEP 2: Compute First-order and Total-order Sobol indices, store YABi

for i = 1:d
    ABi = A;
    ABi(:, i) = B(:, i);
    %YAB(:,i) = model_func(ABi);
    parfor k = 1:N
        YAB(k,i) = model_func(ABi(k,:));
    end
    %first_order(i) = (1/N)*sum(YB .* (YAB(:,i) - YA)) / VY;
    first_order(i) = (mean(YA .* YAB(:,i)) - mean(YA)*mean(YAB(:,i))) / VY;
    total_order(i) = (1/(2*N))*sum((YA - YAB(:,i)).^2) / VY;
end

end
