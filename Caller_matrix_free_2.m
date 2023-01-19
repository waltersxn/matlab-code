clear;
clc;
set(groot,'defaultLineLineWidth',2.0)


L = 3;
M = 30; % 0
hx = L/M;
H = 2;
N = 20; % 0
hy = H/N;
D = -(2/(hx^2) + 2/(hy^2));
Dinv = 1/D;
tol = 0.1;

U0 = zeros(M+1,N+1);
% build and populate b
b = zeros(M+1,N+1);
for i = 1:M+1
    for j = 1:N+1
        b(i,j) = f_gee((i-1)*hx, (j-1)*hy);
    end
end

[U,E] = f_Jacobi(@f_A, U0, b, tol, Dinv);

[Us, Es] = f_steep(@f_A, U0, b, tol);


semilogy(Es)
hold on
semilogy(E)
legend("Steepest descent method", "Jacobi method")
title("Error in matrix-free differential approximation")
hold off

% function to implement the Jacobi method
function [Ucurr, Errs] = f_Jacobi(fA, U0, b, tol, Dinv)
% where Dinv is the inverted diagonal element Dinv = 1/D
iter_max = 500;
Errs = zeros(iter_max, 1);
Ucurr = U0;

for i = 1:iter_max

r = b - fA(Ucurr);
Ucurr = Ucurr + Dinv * r;
Errs(i) = norm(r);
if (norm(r) < tol*norm(b))
    Errs = Errs(1:i);
    return
end
end 
% if you make it through the loop
disp("Jacobi not converged")
end


% function to implement steepest descent

function [Ucurr, Errs] = f_steep(fA, U0, b, tol)

Ucurr = U0;
iter_max = 500;
Errs = zeros(iter_max, 1);

for i = 1:iter_max
    r = b - fA(Ucurr);
    alpha = dot(r, r) / dot(r, fA(r));
    %p = (transpose(r)*r / transpose(r) * fA(r)) * r;
    p = alpha * r;
    Ucurr = Ucurr + p;
    Errs(i) = norm(p);
    if (norm(r) < tol* norm(b))
        Errs = Errs(1:i);
        return
    end
end
disp("Steepest descent not converged")
end



% your g(x,y) function
function u = f_gee(x,y)
H = 2;
L = 3;
u = -6*x*(y^2)*(H-y) + 2*x*(L-x^2)*(H-3*y);
end

% your matrix-free matrix
function Uout = f_A(Uin)

L = 3;
M = 30; % 0
hx = L/M;
H = 2;
N = 20; % 0
hy = H/N;

[m,n] = size(Uin);
Uout = zeros(m,n);

for i = 2:(m-1)
    for j = 2:(n-1)
        Uout(i,j) = ((Uin(i-1,j) + Uin(i+1,j)) / (hx^2)) ...
                   + ((Uin(i,j-1) + Uin(i,j+1)) / (hy^2)) ...
                   - (2/(hx^2) + 2/(hy^2))*Uin(i,j);
    end
end

end




