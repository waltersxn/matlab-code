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


U0 = zeros(N+1,M+1);
% store b as an (M+1)x(N+1) dim array
b = zeros(N+1,M+1);
for i = 1:N+1
    for j = 1:M+1
        b(i,j) = -f_gee((i-1)*hy, (j-1)*hx, L, H);
    end
end
b;
%dyb = f_dy(b,hy)
% identity should return original
%Ib = f_I(b)

%[U, Errs] = f_Jacobi(@f_A, U0, b, tol);
%plot(Errs(1:50))


% activate the Jacobi method
Ucurr = U0;
iter_max = 200; %% CHANGE ITER MAX WHEN READY
Errs = zeros(iter_max,1);
for i = 1:iter_max
    r = f_A(Ucurr, hx, hy) - b;
    Errs(i) = norm(r);
    if (norm(r) < tol*norm(b))
        Errs = Errs(1:i);
        break
    end % otherwise
    % 1/D is equiv to inv(D) where D is a diag matrix
    Ucurr = Ucurr + r*(1/D);
end
if (length(Errs) == iter_max) 
    disp("Jacobi not converged!!")
end



% forcing function to populate the b vector
function q = f_gee(x, y, L, H)
q = 6*x*(y^2)*(H-y) + 2*x*(L^2-x^2)*(H-3*y);
end

% function to act as an identity matrix
function Uout = f_I(Uin)
[m,n] = size(Uin);
Uout = zeros(m,n);
for i = 1:m
    for j = 1:n
        Uout(i,j) = Uin(i,j);
    end
end
end

function Uout = f_dy(Uin, dy)

[m,n] = size(Uin);
Uout = zeros(m,n);
for i = 2:(m-1)
    for j = 2:(n-1)
        i
        j
        uij1 = Uin(i+1,j)
        uij_1 = Uin(i-1,j)
        Uout(i,j) = (Uin(i+1,j)-Uin(i-1,j))/(2*dy)
        uout = Uout(i,j)
    end
end

end

% function to simulate the action of A
function Uout = f_A(Uin, hx, hy)
[n, m] = size(Uin);
% note that m = M+1, n = N+1 by construction of U0
disp("Uin is " + Uin)
Uout = zeros(n, m);
    for i = 2:(n-1)
        for j = 2:(m-1)

            %i
            %j
            %(Uin(i-1,j)/hx^2)
            %(Uin(i+1,j)/hx^2)
            %(Uin(i,j-1)/hy^2)
            %(Uin(i,j+1)/hy^2)
            %- Uin(i,j)*(2/hx^2 + 2/hy^2)

            Uout(i,j) = (Uin(i,j-1)/hx^2) + (Uin(i,j+1)/hx^2) ...
                + (Uin(i-1,j)/hy^2) + (Uin(i+1,j)/hy^2) ...
                - Uin(i,j)*(2/hx^2 + 2/hy^2);
  
            %Uout(i,j)
        end
    end
    disp("Uout is " + Uout)
end



