function dg = f_loan_der(rate)
    amount = 100000;
    payment = 10000;
    n = 20;
    dg = n*(1+rate).^(n-1) - ...
        (payment/amount).*((rate+1).^n.*(n*rate-rate-1)+1)./rate.^2;
    
    der_loan(rate == 0) = -payment;
end