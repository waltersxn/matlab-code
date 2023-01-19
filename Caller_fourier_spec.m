clear;
clc;

% initial arguments given by the problem statement
N = 64;
delta_t = 0.1;
t = 10;

h = 2*pi/N;



X = zeros(N, 1);
for i = 1:N
    X(i) = i*h - pi;
end

Xo = exp(sin(X));

X_ex = exp(sin(X - t));
X_ap = f_fou_spec_update(t, N, delta_t, Xo);

plot(X, X_ex, '.', X, X_ap, 'o')
legend("Exact solution", "Approximated solution")
xlabel("Position [-pi, pi]")
ylabel("U(t_{10})")
title("Fourier spectral differentiation solution of advection equation")
