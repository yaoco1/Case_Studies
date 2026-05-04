function A = Poisson2D(n)

% Grid step size
h = 1/(n+1);

% 1D Laplacian
e = ones(n,1);
T = spdiags([-e 2*e -e], -1:1, n, n);

% Identity Matrix
I = speye(n);

% 2D Laplacian (Kronecker)
A = kron(I, T) + kron(T, I);
A = A / h^2;

end