clear;
clc;

a = -2;
b = 2;
tol = 1e-6;

[z, E] = f_super_rooter2(@poly5, a, b, tol);

Nsr = 1:length(E);

set(groot,'defaultLineLineWidth',2.0)

semilogy(Nsr, E, 'r')
hold on
title("Errors reduction in SuperRoot method")
xlabel('number of iterations')
ylabel('error from final zero')
legend("SuperRoot")
hold off

function f = poly5(x)
  f = x.^5-x+1;

end