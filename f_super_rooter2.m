function [z, errors] = f_super_rooter2(func, x, y, tol)
    
    
    a = x; % left bracket
    b = y; % right bracket for bisection
    xold = a; % x old for secant
    xcur = b; % x current for secant
    
    
    % fa = fold;
    % fb = fcur; % will need these for later
    % bracket = abs(b - a); % length of the first interval
    iter_max = 100; % some reasonable expectation of convergence
    errors = zeros(1,iter_max);
    
    for k = 1:iter_max
        fold = func(xold); % need to do this once
        fcur = func(xcur);
        % calculate the secant update
        
        secf = (fcur - fold) / (xcur - xold);
        sec_step = fcur / secf;
        xtemp = xcur - sec_step; % the first guess of z
        
        % check inside bracket
        if (xtemp >= a && xtemp <= b)
           
            %make a temporary interval
            c = a;
            d = b;
            ftemp = func(xtemp);
            if sign(ftemp) == sign(func(c))
                c = xtemp;
            else
                d = xtemp;
            end
            
            % calculate bracketing factor of new interval
            brack_fac = (b - a) / (d - c);
            if (brack_fac >= 2) % use the secant update
                z = xtemp;
                % update your bracket
                a = c;
                b = d;
                % print current interval and guess
                fprintf('%6i %20.10f %20.10f %20.10f %20.10s\n', ...
                k,xold,xtemp,xcur, "secant")
                % store the error bracket
                errors(k) = abs(b - a);
                if (abs(sec_step) < tol)
                    errors = errors(1:k);
                    return
                end % you're done 
                % otherwise, update
                xold = xcur;
                xcur = xtemp;
            else
                p = (a + b) / 2;
                fp = func(p);
                z = p;
                % display interval and current guess
                fprintf('%6i %20.10f %20.10f %20.10f %20.10s\n', ...
                k,a,z,b, "bisection")
                % update your bracket
                if (sign(func(a)) == sign(fp))
                    a = p;
                else
                    b = p;
                end
                % check if within tolerance
                bracket = abs(b-a);
                errors(k) = bracket;
                if (bracket < tol)
                    errors = errors(1:k);
                    return
                end % you're done 
                % otherwise, update the secant x's
                xold = xcur;
                xcur = p;
            end % if bracketing factor

        else % if xtemp is NOT inside bracket, bisect
            p = (a + b) / 2;
            fp = func(p);
            z = p;
            % display interval and current guess
                fprintf('%6i %20.10f %20.10f %20.10f %20.10s\n', ...
                k,a,z,b, "bisection")
            % update your bracket
            if (sign(func(a)) == sign(fp))
                a = p;
            else
                b = p;
            end
            % check if within tolerance
            bracket = abs(b-a);
            errors(k) = bracket;
            if (bracket < tol)
                errors = errors(1:k);
                return
            end % you're done
            % otherwise, update the secant x's
            xold = xcur;
            xcur = p;
        end
        
    end % for loop

end % super rooter