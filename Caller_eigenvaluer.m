clear;
clc;

N = 64;
e = ones(N-1,1);
h = pi/N;
a = [e -2*e e]/h^2;
A = spdiags(a, [-1, 0 1], N-1, N-1);

tol = 1e-6;

x = h*[1:N-1]';

% case m == 1
Y1 = zeros(N,1);
for i = 1:N-1
    Y1(i) = i*h;
end
Y1(N) = -1;

% case m == 2
Y2 = zeros(N,1);
for i = 1:N-1
    Y2(i) = 2*i*h;
end
Y2(N) = -4;

% case m == 3
Y3 = zeros(N,1);
for i = 1:N-1
    Y3(i) = 3*i*h;
end
Y3(N) = -9;

% case m == 4
Y4 = zeros(N,1);
for i = 1:N-1
    Y4(i) = 4*i*h;
end
Y4(N) = -16;



z1 = f_newton_vec(@f_eigenvaluer, @f_eigen_jac, Y1, tol);
z2 = f_newton_vec(@f_eigenvaluer, @f_eigen_jac, Y2, tol);
z3 = f_newton_vec(@f_eigenvaluer, @f_eigen_jac, Y3, tol);
z4 = f_newton_vec(@f_eigenvaluer, @f_eigen_jac, Y4, tol);


lambda1 = z1(N);
evec1 = z1(1:N-1);

lambda2 = z2(N);
evec2 = z2(1:N-1);

lambda3 = z3(N);
evec3 = z3(1:N-1);

lambda4 = z4(N);
evec4 = z4(1:N-1);

plot(x, evec1, '.', x, evec2, 'o', x, evec3, 'x', x, evec4, '*')
legend("lambda = -0.99", "lambda = -3.99", "lambda = -8.98", "lambda = -15.95")
title("Plotting eigenvectors of a tridiagonal matrix")

