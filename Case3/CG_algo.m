%Q1

function [x, iter, rel_res] = CG_algo(A, b, x0, tol, maxIter)

% Initialization
x = x0;
r = b - A*x;
p = r;
normb = norm(b);
rel_res = [];
res_norm = norm(r) / normb;
rel_res = [rel_res; res_norm];
iter = 0;

% CG iterations
while res_norm > tol && iter < maxIter
    Ap = A * p;
    alpha = (r' * r) / (p' * Ap);
    x = x + alpha * p;
    r_new = r - alpha * Ap;
    
    % Calculate the new residuals
    res_norm = norm(r_new) / normb;
    rel_res = [rel_res; res_norm];
    
    % Covergence
    if res_norm < tol
        break;
    end
    
    beta = (r_new' * r_new) / (r' * r);
    
    p = r_new + beta * p;
    
    r = r_new;
    iter = iter + 1;
end

end