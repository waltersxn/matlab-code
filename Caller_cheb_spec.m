clear;
clc;
set(groot,'defaultLineLineWidth',2.0)

tol = 1e-9;
a = -10;
% first, your two diff matrices and their spaced grids

[Ds, xs] = f_cheb(16); 
[Db, xb] = f_cheb(512);

% initial guesses for the two sizes
U0s = f_u0(16);
U0b = f_u0(512);

% solution vectors for the two sizes
Ys = f_newton_vec2(@f_cheb_spec, U0s, tol);
Yb = f_newton_vec2(@f_cheb_spec, U0b, tol);

plot (xs, Ys, xb, Yb);
legend("N = 16", "N = 512")
title("Using Chebyschev Spectral Matrix to solve diff eqs")


% plot the errors for 2^k points
% using N = 512 (Ub) as the "exact" solution

K = [1,2,3,4,5,6];
Err = zeros(length(K), 1);

for i = 1: length(K)
    N = 2^K(i); % this is this N
    h = 512 / N; % step size difference between this N and 512

    U0 = f_u0(N); % the initial vector for the solution in question
    Ut = f_newton_vec2(@f_cheb_spec, U0, tol);
    sum = 0;
    for j = 2:N+1 % the sol vector will have this length
        
        % sum the squares of differences of common indexed entries past 1
        sum = sum + (Ut(j) - Yb( (j-1)*h))^2;
    end
    Err(i) = sqrt((sum + (Ut(1)- Yb(1))^2) / (N+1));
end

semilogy(K,Err);
xlabel("Number of points in 2^k +1 ")
ylabel("Error")
title("Plotting spectral accuracy of the Chebyschev "...
    + "spectral differentiation method")


% returns an initial guess vector for a desired N
% note that the vector will have length N+1
function U0 = f_u0(N)
[~,x] = f_cheb(N);
U0 = (1-x.^2) /2;

end

function [Y, YJac] = f_cheb_spec(X)
    
a = -10;

    n = length(X); % the # of intervals
    [D,~] = f_cheb(n-1);
    
    Y = zeros(n,1);
    % Y(1) = u_0 is left as 0
    for j = 1:n
        Y(n) = Y(n) + D(n,j) * X(j); %populate u_n
    end
    Y(n) = Y(n) + 1;
    
    D2 = D*D; % the second der matrix
    for i = 2:(n-1) % only for interior points
        for j = 1:n
            Y(i) = Y(i) + D2(i,j)*X(j);
        end
        Y(i) = Y(i) + exp(a*X(i));
    end
    YJac = D2;
    for j = 2:n-1
        % populate Jacobian with the diag{...}
        YJac(j,j) = YJac(j,j) + a*exp(a*X(j));     
    end
    for j = 1:n
        % set the top row to be all zeros
        YJac(1,j) = 0;
        % bottom row to be the first der matrix bottom row
        YJac(n,j) = D(n,j);
    end
    % top left corner is 1
    YJac(1,1) = 1;
    
end