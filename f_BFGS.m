function xcurr = f_BFGS(phi, grad_phi, x0, tol)

xcurr = x0;
% starting with the identity matrix
Gcurr = eye(length(xcurr));
% take the gradient value at first guess of x
dfcurr = grad_phi(xcurr);

I = eye(length(xcurr));

for k = 1:100
    k
    p = -Gcurr * dfcurr;
    alpha = fminbnd(@(a)phi(xcurr - a*dfcurr), 0,5);
           
    w = alpha * p;

    if (norm(w) < tol)
        disp("BFGS converged in " + k + " steps.")
        return
    end
    % otherwise, update x
    xcurr = xcurr + w;
    dfnew = grad_phi(xcurr);
    y = dfnew - dfcurr

    Gcurr = (I - w*transpose(y) / (transpose(y)*w))*Gcurr...
        *(I - y*transpose(w)/(transpose(y)*w))...
        + w*transpose(w)/(transpose(y)*w)

    dfcurr = dfnew;


end



end