function out = f_shooterg(alpha)

% function g for shooting problem
% given some alpha, solves for y1 and y2
% returns a value that attempts to satisfy 
% the continuity condition from the problem statement
dt = 1/100;
N = 1/dt;

k = 3;
ycurr = [1;alpha];
for n = 1:N

    D = Dmat(n*dt, k);
    K1 = -D*ycurr;
    K2 = -D*(ycurr + 0.5 * dt *K1);
    K3 = -D*(ycurr + 0.5 * dt *K2);
    K4 = -D*(ycurr +       dt *K3);

    ycurr = ycurr + dt/6 * (K1 + 2*K2 + 2*K3 + K4);
    %ycurr = ycurr + dt*K1;
end
out = ycurr(2)/ycurr(1) + k;

% function for D matrix at each iteration
function D = Dmat(tn, a)

D = zeros(2);
D(1,2) = 1/f_m(tn);
D(2,1) = a^2*f_m(tn);
end

% function m given by the problem
function m = f_m(x)
    if x <= 0
        disp("t cannot be <= 0")
        m = nan;
        return
    elseif 0<x && x<1
        m = 0.5 + 0.4*sin(12*pi*x);
        return
    else
        %disp("x == 1")
        m = 1;
    end
end
end