function f = f_eigenvaluer(Y)

% populates a size-65 output vector for the eigenvalue solution function
% f(y,lambda) where y is a desired eigenvector of length 64 and lambda 
% is the desired eigenvalue

N = 64;
e = ones(N-1,1);
h = pi/N;
a = [e -2*e e]/h^2;
A = spdiags(a, [-1, 0 1], N-1, N-1);

if (length(Y) ~= 64)
    disp("Input vector must have lenght N+1 = 64!!");
    return;
end

f = zeros(N,1);

% populate the first N entries, these are the Ay - lambday of f(y,lambda)
for i = 1:N-1
    sum = 0;
    for k = 1:N-1
       sum = sum + A(i,k)*Y(k); 
    end
    f(i) = sum - Y(N)*Y(i);
end

% final entry, accounting for the y^Ty = 1 restriction
sum2 = 0;
for i = 1:N-1
    sum2 = sum2 + Y(i)^2;
end
f(N) = sum2 - 1;

end