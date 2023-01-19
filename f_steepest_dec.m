function xcurr = f_steepest_dec(phi, grad_phi, x0, tol)
% finds the minimum of a multivariable function,

% phi:: the multivariable function in question
% :::: evaluated discreetly on internal points only
% grad_phi: the gradient of the function in question
%   Detailed explanation goes here


iter_max = 100;

xcurr = x0; % your initial guess

for k = 1:iter_max

    dfcurr = grad_phi(xcurr);

    alpha = fminbnd(@(a)phi(xcurr - a*dfcurr),0,5);

    p = - (alpha * dfcurr);
    xcurr = (xcurr + p);

    if (norm(p) < tol)
        disp("Steepest converged in " + k + " steps.")
        return
    end


end

disp("Steepest descent not converged!!");

end