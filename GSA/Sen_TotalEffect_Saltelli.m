function sen_vector = Sen_TotalEffect_Saltelli( A, B, f )
% This function compute the total effects Sobol' index by the single loop
% method proposed by Saltelli [1]
% n: number of sample
% k: number of dimension
% A: n*k matrix, quasi-random sample recommended
% B: n*k matrix, quasi-random sample recommended
% you can generate quasi-random sample using Latin Hypercube sampling 
% f: function handle, so y_A = f(A), y_B = f(B)
% Note: the computational cost of this functio is 2*n + n*k
% Ref: [1] Saltelli, A., Ratto, M., Andres, et. al, "Global sensitivity 
% analysis: the primer", 2008, page 164

y_A = f(A); y_A = y_A(:);
y_B = f(B); y_B = y_B(:);
[N, k] = size(A);

yAll = [y_A; y_B];
f_0  = mean(yAll);
VarY = mean(yAll.^2) - f_0^2;

sen_vector = zeros(1, k);

for j = 1:k
    % Hybrid matrix A_Bj: A with column j replaced by B's column j
    A_Bj = A;
    A_Bj(:, j) = B(:, j);
    y_ABj = f(A_Bj); y_ABj = y_ABj(:);

    % Total-effect (Jansen 1999): pair yA with y(A_Bj)
    sen_vector(j) = mean((y_A - y_ABj).^2) / (2 * VarY);

    fprintf('%d out of %d indices finished: %.6f\n', j, k, sen_vector(j));
end
   
end