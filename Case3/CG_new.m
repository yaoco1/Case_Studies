% Q4

% Initialization
n = 1024;
x0 = zeros(n,1);
tol = 1e-8;
maxIter = 5000;

[A, b] = get_covariance_matrix(n, 100, 0.05);

[x1, iter1, rel_res1] = CG_algo(A, b, x0, tol, maxIter);

M = chol(A);
A_new = M\(A/M');

% New CG
[x2, iter2, rel_res2] = CG_algo(A_new, b, x0, tol, maxIter);

fprintf('Original iterations: %d\n', iter1);
fprintf('Preconditioned iterations: %d\n', iter2);

% Plotting
figure;
semilogy(rel_res1, '-o'); 
hold on;
semilogy(rel_res2, '-x');
legend('Original A', 'Preconditioned A_{new}');
xlabel('Iteration');
ylabel('Relative Residual');
title('CG Convergence Comparison');
grid on;

eig_A = eig(A);
eig_Anew = eig(A_new);

figure;
plot(sort(eig_A), 'o');
hold on;
plot(sort(eig_Anew), 'x');
legend('A', 'A_{new}');
title('Eigenvalue Comparison');
xlabel('Index');
ylabel('Eigenvalue');
grid on;