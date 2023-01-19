
% returns the estimate after n iterations
% and the relative error wrt the actual sine function
function [f,e] = f_taylor_sin(x,n)
    
    f = 0;
    for k = 1:n
        f = f + ((-1)^(k-1))*(x^(2*k-1)) / factorial(2*k-1);
    end

    e = (abs((f - sin(x))/sin(x)))*100;
end