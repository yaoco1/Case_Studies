clc; clear;

% Parameters
n = 31;                 % 2^p - 1
direct_N = 16;          
weight = 1;             % scaling
pre_smooth = 5;
pos_smooth = 5;

% Construct a system
A = Poisson2D(n);
f = randn(n^2,1);

% Construct multigrid hierarchy
[A_list, R_list, max_level] = Restrict(A, direct_N);

% Initial solution
x = zeros(n^2,1);

% V-cycle
num_cycles = 13;
residuals = zeros(num_cycles,1);

for k = 1:num_cycles
    x = Vcycle(1, A_list, R_list, f, x, direct_N, weight, pre_smooth, pos_smooth);

    r = f - A*x;
    residuals(k) = norm(r,2);
end

% Draw a diagram
figure;
semilogy(1:num_cycles, residuals, '-o','LineWidth',2);
xlabel('V-cycle iteration k');
ylabel('||r||_2');
yticks([0.1 0.2 0.5 1 2 5])
yticklabels({'0.1','0.2','0.5','1','2','5'})
title('Multigrid V-cycle Convergence');
grid on;