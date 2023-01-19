function y = f_heat_eq_exact(x)
%f_heat_eq_exat 
% the closed-form calculated function for the heat equation 
% in problem 1 of HW2
    if (x < 0) | (x > pi) 
        disp("The function is only defined for values of x" + ...
            "between 0 and pi. Attempted: ")
        disp(x);
        y = NaN;
    else
        y = sin(x) + x./pi;
    end

end