clear;
clc;
set(groot,'defaultLineLineWidth',2.0)

% first, set all your values
L = 3;
H = 2;
M = 30;
N = 20;
hx = L/M;
hy = H/N;
D = -(2/hx^2 + 2/hy^2); % the diagonal value of conceptual A
tol = 0.1;


U0 = zeros(M+1,N+1);
% store b as an (M+1)x(N+1) dim array
b = zeros(M+1,N+1);
for i = 1:M+1
    for j = 1:N+1
        b(i,j) = -f_gee((i-1)*hx, (j-1)*hy, L, H);
    end
end
b

[U, Errs] = f_Jacobi(@f_A, U0, b, tol);
%plot(Errs(1:50))


% function to action Jacobi method 
function [Ucurr, E] = f_Jacobi(fA, U0, b, tol)
L = 3;
H = 2;
M = 30;
N = 20;
hx = L/M;
hy = H/N;
D = -(2/hx^2 + 2/hy^2); % the diagonal value of conceptual A
Ucurr = U0;
iter_max = 1000;
E = zeros(iter_max,1);
for i = 1:iter_max
    r = b - fA(Ucurr);
    E(i) = norm(r);
    if (norm(r) < tol*norm(b))
        E = E(1:i);
        return
    end % otherwise
    % 1/D is equiv to inv(D) where D is a diag matrix
    Ucurr = Ucurr + r*(1/D);    
end
disp("Jacobi not converged!!")

end



% forcing function to populate the b vector
function q = f_gee(x, y, L, H)
q = 6*x*(y^2)*(H-y) + 2*x*(L^2-x^2)*(H-3*y);
end

% function to simulate the action of A
function Uout = f_A(Uin) 
L = 3;
H = 2;
M = 30;
N = 20;
hx = L/M;
hy = H/N;

% note that m = M+1, n = N+1 by construction of U0
Uout = zeros(M+1, N+1);
    for i = 2:M
        for j = 2:N
            Uout(i,j) = (Uin(i-1,j)/hx^2) + (Uin(i+1,j)/hx^2) ...
                + (Uin(i,j-1)/hy^2) + (Uin(i,j+1)/hy^2) ...
                - Uin(i,j)*(2/hx^2 + 2/hy^2);
        end
    end
end
