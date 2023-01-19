function [x, errors] = f_newton_zero(func, dfunc, x0, tol)
    iter_max = 100;
    x = x0;
    errors = (1:iter_max);
    % make iter_max for size of array and # of iterations

    for k = 1:iter_max
        fx = func(x);
        dfdx = dfunc(x);
        dx = -fx/dfdx;
        x = x+dx;
        errors(k) = abs(dx);
        fprintf('%6i %20.10f %20.10f %20.10f\n', k,x,dx,fx) 
        if abs(dx)<tol
            errors = errors(1:k); % chop errors down to only what you need
            return
        end
    end

    disp('Newton not converged, possibly incorrect solution.')
end