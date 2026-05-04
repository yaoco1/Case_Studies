function x = Vcycle(level, A_list, R_list, b, x0, direct_n, weight, pre_smooth, pos_smooth)

% level      Current level (initial=1)
% A_list     Array of coefficient matrices on each level
% R_list     Array of restriction operators on each level
% b          Right-hand side
% x0         Initial Guess
% direct_n   Threshold for directly solving Ax = b
% weight     Scaling coefficient between restriction and prolongation
%            operators
% pre_steps  Number of iters of pre-smoothing
% pos_steps  Number of iters of post-smoothing


	
	% Load coefficient matrix
	A = cell2mat(A_list(level));
	
	%Solve directly if problem is small enough
	n =size(b, 1);
	if (n <= direct_n)
		x = A \ b;
		return;
	end

	% Pre-smoothing
	%x = 'your weighted Jacobi smoother'; ** input needed **
	x = x0;
    D = diag(diag(A));
    for k = 1:pre_smooth
        x = x + (2/3) * (D \ (b - A*x));
    end

	% Load restriction operator and construct interpolation operator
	R = cell2mat(R_list(level));
	P =  R' * weight;       % makes sure the coarse grid has the correct scaling when interpolated back to the fine grid. 
	coarse_n = size(R, 1);
	
	% Compute residual and transfer to coarse grid
	% ** input needed **
    r = b - A*x;      % residual on fine grid
    r_H = R * r;      % restrict to coarse grid
	
	% Implement the V-cycle recursively
	x0  = zeros(coarse_n, 1);
	e = Vcycle(level + 1, A_list, R_list, r_H, x0, direct_n, weight, pre_smooth, pos_smooth);
	
	% Transfer error to fine grid and correct
	% ** input needed **
	x = x + P * e;

	% Post-smoothing
	%x = 'your weighted Jacobi smoother'; ** input needed **
    for k = 1:pos_smooth
        x = x + (2/3) * (D \ (b - A*x));
    end
    
end