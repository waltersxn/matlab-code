% bracket values for bisection method
bi1 = 0.07; 
bi2 = 0.08;

%starting value for newton method
xn = 0.08;

% starting values for secant method
xc = 0.08;
xo = 0.1;
tol = 1e-6;

[zbi, errorsbi] = f_bisection_zero(@f_loan, bi1, bi2, tol);
[znw, errorsnw] = f_newton_zero(@f_loan, @f_loan_der, xn, tol);
[zsc, errorssc] = f_secant_zero(@f_loan, xc, xo, tol);

%errorsnw = [errorsnw 0];
%errorssc = [errorssc 0];

Nbi = 1:length(errorsbi);
Nnw = 1:length(errorsnw);
Nsc = 1:length(errorssc);

set(groot,'defaultLineLineWidth',2.0)

semilogy(Nbi, errorsbi, 'g')
hold on
semilogy(Nnw, errorsnw, 'rx')
semilogy(Nsc, errorssc, 'bo')
title("Errors reduction in zero-finding methods")
xlabel('number of iterations')
ylabel('error from final zero')
legend("Bisection", "Newton", "Secant")
hold off