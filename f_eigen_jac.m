function M = f_eigen_jac(Y)


% using the same matrix as in the original function
N = 64;
e = ones(N-1,1);
h = pi/N;
a = [e -2*e e]/h^2;
A = spdiags(a, [-1, 0 1], N-1, N-1);

if (length(Y) ~= N)
    disp("Input vector must have lenght N!!");
    return;
end

lambda = Y(N);
M = A;
% populate the Jacobian matrix of the system
for i = 1:N-1
   M(i,i) = M(i,i) - lambda; % the diagonal is all aii - lambda
   M(i,N) = -Y(i); % the Nth column
   M(N,i) = 2*Y(i); % the bottom row
   M(N,N) = 0; % the bottom-right corner
end

end