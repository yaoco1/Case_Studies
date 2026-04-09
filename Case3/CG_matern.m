% Q2

% Initialization
n = 1024;
tau = 20;    % The hyperparameter
noise = 0.005;
x0 = zeros(n,1);
tol = 1e-8;   % Tolerance 
maxIter = 5000;

[A, b] = get_covariance_matrix(n, tau, noise);

[x, iter, relres] = CG_algo(A, b, x0, tol, maxIter);

fprintf('Iterations: %d\n', iter);

% Plotting
figure;
semilogy(relres, '-o');
xlabel('Iteration');
ylabel('Relative Residual');
title('CG Convergence (n=1024)');
grid on;