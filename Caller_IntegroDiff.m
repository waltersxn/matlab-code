clear;
clc;
set(groot,'defaultLineLineWidth',2.0)

% set problem specs

a = 0.1;
u0 = 1;
dt = 0.01;
T = 0:dt:01;
N = length(T) - 1; % # of intervals or time steps
u = zeros(N+1,1);
u(1) = u0;

% implement algorithm
C = 1 + dt*a - (dt^2)*gee(0)/2;

for k = 1:N

    sum = dt^2*gee(T(k+1))/2 + u(k);
    for j = 1:k
        sum = sum + (dt^2)*gee(T(k+1)-T(j))*u(j);
    end
    u(k+1) = sum / C;
end

plot(T, u)
xlabel("t = Time in seconds")
ylabel("u(t)")
title("Integro-differential equation")


% the g(t) given by the problem
function y = gee(t)
    y = sin(8*pi.*t);
end