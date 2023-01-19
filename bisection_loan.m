x1 = 0.07;
x2 = 0.08;
bracket = 10^(-6);

% calls the bisector function to get a zero of money
[rate, Errors] = bisector(@money_owed, x1, x2, bracket);

rxt = 1:length(Errors);
plot(rxt, Errors, 'r')
%{ 
 function: money_owed_con -- given an annual interest rate, returns 
           the balance owed on a loan after n years of a fixed payment size
 input: rate -- the interest rate, must be >= 0
 note: no other inputs, but other variables can be adjusted internally
 these are: amount-size of initial loan; payment-size of annual payments;
            n-numbers of years/iterations of payment
%}
function q = money_owed(rate)
    amount = 100000;
    payment = 10000;
    n= 20;
    q = (1 + rate)^n - (payment/amount)*((1+rate)^n-1)/rate;
end

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
function [x, err] = bisector(func, a, b, tol)

    fa = func(a);
    fb = func(b);
    
    %first check that inputs are valid

    if sign(fa) == sign(fb) || isnan(fa) || isnan(fb)
        disp('error: root not bracketed or function is NaN at endpoint');
        x = NaN; err = NaN;
        return
    end

    %if the inputs are good, begin bisection iterations
    n = 0; err = [];
    while (b-a > tol)
        n = n+1; %increment the iteration number
        err = [err b-a]; %log the current error interval length
        p = (a+b)/2; %take the midpoint
        x = p;
        fp = func(p); %func value at midpoint
        %formatted printing to watch the bisection converge
        fprintf('%6i %20.10f %20.10f %20.10f\n', n, a, p, b)
        %update the endpoint with same sign as the midpoint
        if sign(fa) == sign(fp)
            a = p;
        else
            b = p;
        end

    end

end






