function D = f_four_spec_diff_mat(N)
% from p2.m, https://people.maths.ox.ac.uk/trefethen/spectral.html
% Lloyd N. Trefethen, Spectral Methods in MATLAB, SIAM, Philadelphia, 2000
% N = number of grid points

h = 2*pi/N; % grid spacing

column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
D = toeplitz(column,column([1 N:-1:2]));
end

% Use example
% D = f_four_spec_diff_mat(N);
% example of using D to differentiate a function
% u(x) = exp(sin(x));
% u = exp(sin(x)); uprime = cos(x).*u;
% plot(x,D*u,'o',x,uprime,'.')
% xlim([-pi pi])