clear;
clc;
set(groot,'defaultLineLineWidth',2.0)

N = 16;

% build your evenly spaced X
Xe = linspace(0,1,N+1);
Xu = sqrt(Xe);

Ie = f_comp_simp(@f_simplef, Xe);
Iu = f_comp_simp(@f_simplef, Xu);

function Y = f_simplef(X)
Y = sin(4*pi*X.^2);
end

% the function that evaluates composite simpsonian integral
% over a given vectorized interval
function I = f_comp_simp(func, X)
% assume X is a vector of points that partition a continuous interval

% your vector of solutions

Y = func(X);

N = length(X) -1;
if (mod(N,2) ~= 0)
    disp("X must have odd length to give an even number of intervals")
    I = nan;
    return
end
% your vector of interval lengths
H = zeros(N, 1);
for j = 2:N+1
    H(j-1) = X(j) - X(j-1);
end

sum_int = 0;
for k = 2:2:N
    sum_int = sum_int + Y(k-1)/6 * (2*H(k-1) +H(k)- (H(k))^2/H(k-1))...
        + Y(k) * (H(k-1)+H(k))^3 / (6*H(k-1)*H(k)) ...
        + Y(k+1) * ( (H(k-1)+H(k))^2/3 - H(k-1)*(H(k-1)+H(k))/2 ) / H(k);
end
I = sum_int;
end
