function [x, errors] = f_secant_zero(func, a, b, tol)

    iter_max = 100;
    x = a; % assigns first value for x
    xo = b; % the previous point
    fxo = func(xo); % value of the function at previous point
    errors = zeros(1,iter_max); % an array to track the error
    
    for k = 1:iter_max
        fxc = func(x); % assigns first value for f(x)

        sec_fx = (fxc - fxo) / (x - xo); % calculate secant deriv. approx.
        dx = fxc/sec_fx; % calc the secant step
        errors(k) = abs(dx); % save this for later
        fprintf('%6i %20.10f %20.10f %20.10f\n', k, x, xo, dx)
        xo = x; % this prepares for the next loop
        x = x - dx; % find the next guess for x
        
        % if you're within the desired tolerance interval
        if abs(dx) < tol
            % take only the submatrix that contains calculated errors
            errors = errors(1:k); 
  
            return;
        end % if

        % otherwise, prepare for the next loop
        fxo = fxc; % the current value of f becomes the old one

    end % for loop
      
end % secant function