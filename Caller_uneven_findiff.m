clear;
clc;
set(groot,'defaultLineLineWidth',2.0)

eps = 0.01;
N = 32;

% the vector of even spacings
Xe = zeros(N+1,1);
% the vector of uneven specings
Xu = zeros(N+1,1);
% the vector of "results"
b = zeros(N+1,1);
b(N+1) = 1;

for i = 1:length(Xe)
    Xe(i) = (i-1)/N;
    Xu(i) = eps*(i-1)/N + (1-eps)*((i-1)/N)^4;
end

% the vector of evenly spaced xi's


Ae = u_mat(Xe);
Au = u_mat(Xu);

Ye = Ae \ b;
Yu = Au \ b;


Xex = linspace(0,1,1000);
Yexp = u(Xex);

subplot(2,1,1)
plot(Xe, Ye, Xu, Yu, Xex, u(Xex))
legend("Even spacing", "Uneven spacing", "Explicit")

% set the zoomed vectors
Z = 0.05;
%evenly spaced
Xez = zeros(N+1,1);
%unevenly spaced
Xuz = zeros(N+1,1);

for i = 1:(N+1)
    Xez(i) = Z*(i-1)/N;
    Xuz(i) = Z*(eps*(i-1)/N + (1-eps)*((i-1)/N)^4);
end

Aez = u_mat(Xez);
Auz = u_mat(Xuz);

Yez = Aez \ b;
Yuz = Auz \ b;

Xexz = linspace(0,0.05,1000);
Yexpz = u(Xexz);
subplot(2,1,2)
plot(Xez, Yez, Xuz, Yuz, Xexz, u(Xexz))
legend("Even spacing", "Uneven spacing", "Explicit")



% given the vector X, constructs the matrix to find Y
function A = u_mat(X)
eps = 0.01;

% the vector of spacings
% note that h(1)=0 and is ignored
h = zeros(length(X), 1);
for i = 2:length(h)
    h(i) = X(i) - X(i-1);
end


% set A and its end diagonal points
A = zeros(length(X),length(X));
A(1,1) = 1;
A(length(X),length(X)) = 1;

for i = 2:(length(X) -1)
    A(i,i-1) = (eps  - h(i+1)) / ((h(i))^2 + h(i)*h(i+1));
    A(i,i) = - (eps + h(i) - h(i+1)) / (h(i)*h(i+1));
    A(i,i+1) = (eps + h(i)) / (h(i)*h(i+1) + (h(i+1))^2);
end

end



% the explicit function u(x)
function y = u(x)
eps = 0.01;
c1 = 1/ (exp(-1/eps) - 1);
y = zeros(length(x),1);

for i = 1:length(x)
    y(i) = c1*exp(-x(i)/eps) -c1;
end
%y = c1*exp(-x/eps) -c1;
end