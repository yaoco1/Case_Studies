% Q3

% Initialization
n = 1024;
tau = 564;    % The hyper-parameter is my student number
tol = 1e-8;   % Tolerance 
maxIter = 5000;
noise_values = logspace(log10(0.5), log10(0.5e-8), 8);      % Noise
iterations = zeros(size(noise_values));
condition_numbers = zeros(size(noise_values));

for i = 1:length(noise_values)
    noise = noise_values(i);
    [A, b] = get_covariance_matrix(n, tau, noise);
    x0 = zeros(n,1);
    [x, iter, rel_res] = CG_algo(A, b, x0, tol, maxIter);
    iterations(i) = iter;
    condition_numbers(i) = condest(A);
    
    fprintf('Noise = %.2e, iteration = %d, condition number = %.2e\n', noise, iter, condition_numbers(i));
end

% Plotting
figure;
semilogx(noise_values, iterations, '-o');
xlabel('Noise');
ylabel('Iterations');
title('Iterations vs Noise');
grid on;

figure;
semilogx(noise_values, condition_numbers, '-o');
xlabel('Noise');
ylabel('Condition Number');
title('Condition Number vs Noise');
grid on;