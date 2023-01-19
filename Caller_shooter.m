clear;
clc;
set(groot,'defaultLineLineWidth',2.0)

% function g for shooting problem
% given some alpha, solves for y1 and y2
% returns a value that attempts to satisfy 
% the continuity condition from the problem statement

% actual function


tol = 1e-6;
X = linspace(1.35, 1.45, 100);
gX = zeros(length(X),1);
for i = 1:length(X)
    gX(i) = f_shooterg(X(i));
end

plot(X, gX);
xlabel("alpha")
ylabel("f_shooterg(alpha)")
title("Plotting the g-function")

[sol, errs] = f_secant_zero(@f_shooterg, 1.38, 1.382, tol);

Q = linspace(0,2,201);
Y = f_shooterY(sol);
plot(Q,Y);
xlabel("x")
ylabel("u(x)")
title("Plotting u(x) on the interval (0,2)")