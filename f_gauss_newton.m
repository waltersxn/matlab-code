function Xcurr = f_gauss_newton(g, Aj, X0, b, tol)
%f_gauss_newton uses the Gauss-Newton method to minimize the residual
% of some nonlinear least squares problem ||g(x) - b||

% g is the function
% A is the dependent matrix derived as dgi/dxk for i = 1,..,m; k = 1,..,n;
% b is the "output" vector
% Xo is the initial guess
% tol is the tolearance

Acurr = Aj(X0);
p = transpose(Acurr)*Acurr \ transpose(Acurr)*b;
Xcurr = X0 + p;
k = 0;

while (norm(p) >= tol*(norm(Xcurr) + 1))
    
    Acurr = Aj(Xcurr);
    bg = b - g(Xcurr);
    p = transpose(Acurr)*Acurr \ transpose(Acurr)*bg;
    Xcurr = Xcurr + p; % this might be a minus!!!
    k = k  +1;

end

disp("Gauss-Newton converged in " + k + " steps.")
% and whatever Xcurr settles on is what you want!!

end