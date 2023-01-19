%{
given a fucntion func, the bisector method solves for some root
func(x) = 0 over an interval (a,b) to within a given error interval
input: -funtion func 
       -endpoints a,b
       -allowed error tol
output: -root of function x
        -error of result err
note: func(a) and func(b) must be valid numbers and have opposite signs
%}
function [x, errors] = f_bisection_zero(func, a, b, tol)

    fa = func(a);
    fb = func(b);
    
    %first check that inputs are valid

    if sign(fa) == sign(fb) || isnan(fa) || isnan(fb)
        disp('error: root not bracketed or function is NaN at endpoint');
        x = NaN; errors = NaN;
        return
    end

    %if the inputs are good, begin bisection iterations
    errors = b - a ; % the first error is just the interval length
    iter_count = 0;
    while (b-a > tol)
        iter_count = iter_count + 1; %increment the iteration number
        
        p = (a+b)/2; %take the midpoint
        x = p;
        fp = func(p); %func value at midpoint
        %formatted printing to watch the bisection converge
        fprintf('%6i %20.10f %20.10f %20.10f\n', iter_count, a, p, b)
        %update the endpoint with same sign as the midpoint
        if sign(fa) == sign(fp)
            a = p;
        else
            b = p;
        end
        errors = [errors b-a]; %log the current error interval length

    end

end