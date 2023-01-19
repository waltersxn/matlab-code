function x = f_super_rooter(func, xoa, xcb, tol)

    a = xoa;
    b = xcb;
    xo = xoa;
    xc = xcb;

    foa = func(xoa)
    bracket = abs(xcb - xoa);
    k = 1; % use for tracking iterations
    while (bracket > tol)
        k % printing just for debugging
        fcb = func(xcb)
        % secant = (f(b) - f(a))/(b-a); take neg of that
        
        sec = (fcb - foa) / (xcb - xoa);
        sec_step = fcb / sec;
        xtemp = xcb - sec_step
        
        % check the update in the current interval
        if (xtemp < xcb && xtemp > xoa)
            % if so, compute the bracket factor
            brac_fac = abs(xcb - xoa) / abs(xtemp - xcb);
            
            if (brac_fac > 2) % so better than bisection
                x = xtemp;
                % print current interval and guess
                fprintf('%6i %20.10f %20.10f %20.10f %20.10s\n', ...
                k,xoa,xtemp,xcb, "secant")
                % use function VALUES!!!!
                if sign(xcb) == sign(xtemp)
                    xcb = xtemp;
                else
                    xoa = xtemp;
                end
               
                xoa = xcb; % current x becomes old x
                xcb = xtemp; % temp x becomes current x
                bracket = abs(xcb-xoa)
                % if the bracket is small enough, we're done
                if (bracket < tol)
                    return
                end

            end % brac_fac if
                
        else       
            % take the midpoint
            x = (xcb + xoa) / 2;
            fprintf('%6i %20.10f %20.10f %20.10f %20.10s\n', ...
            k,xoa,x,xcb, "bisection")
            bracket = abs(xcb - x)
            if (bracket < tol)
                return
            end % we're done
            % otherwise, do bisection comparison
            fp = func(x);
            if sign(fp) == sign(foa)
                xoa = x;
            else 
                xcb = x;
            end % sign if

        end % bracket if
    foa = fcb  
    k = k+1;
    
    end % while loop
end % function