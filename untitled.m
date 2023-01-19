clear;
clc;

[z1, e1] = f_heat_eq_es(512);

M = length(z1);
eqN = 1:M;
X = zeros(M, 1);
for i = 1:M
    X(i) = i*pi/M;
end

plot(X, z1, 'r')
hold on
plot(X, f_heat_eq_exact(X), 'g')
hold off

function [z, e] = f_heat_eq_es(N)
%f_Heat_Eq Summary of this function goes here
%   Detailed explanation goes here
T0 = 0;
Tpi = 1;

h = pi/N;

% partition the rod, assign values
X = transpose(zeros(1,N-1));
for i = 1:(N-1)
    X(i) = i*pi/N;
end

% using formula that for each interior i,
% -(1/h^2)*T_(i-1) + (2/h^2)*T_i -(1/h^2)*T_(i+1) = sin(x_i)
% where T_i = T(x_i)

% create and populate the matrix A 
A = zeros(N-1,N-1);

%boundary values
A(1,1) = 2/(h^2);
A(1,2) = -1/(h^2);
A(N-1,N-2) = -1/(h^2);
A(N-1,N-1) = 2/(h^2);

% diagonal and paradiagonals
for i = 2:N-2  
    A(i,i) = 2/(h^2);
    A(i,i+1) = -1/(h^2);
    A(i,i-1) = -1/(h^2);
end

% create and populate vector b

b = sin(X);
b(1) = b(1) + T0/(h^2);
b(N-1) = b(N-1) + Tpi/(h^2);

% z is the vector of solutions
z = A\b;

% e is the error of the solution
temp_sum = 0;
for i = 1:N-1
    temp_sum = temp_sum + (z(i) - f_heat_eq_exact(X(i)))^2;
end
e = sqrt(temp_sum/(N-1));


%temp_sm = (z - f_heat_eq_exact(X)).^2;
%temp_sum = sum(temp_sm);
%e = sqrt(temp_sum/(N-1));

end
