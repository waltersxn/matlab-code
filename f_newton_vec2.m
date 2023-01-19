function X = f_newton_vec2(func, X0, tol)
    iter_max = 100;
    X = X0;
    % errors = (1:iter_max);
    % make iter_max for size of array and # of iterations
    % j_func should be a function that returns the jacobian of func
    % tol must be a vector size same as X0 that is all tolerances

    for k = 1:iter_max
        % assign values
        [fx,J] = func(X);

        % calculate the update
        p = J \ -fx;
        X = X + p;

        % errors(k) = abs(dx); ignore errors for this one
        % printf('%6i %20.10f %20.10f \n', k,J,p) 
        if norm(p) < tol * (norm(X) + 1)
            % errors = errors(1:k); % chop errors down to only what you need
            disp("Newton converged in " + k + " steps")
            return
        end
    end
    
    % if no solution
    disp('Newton not converged, possibly incorrect solution.')
end