clear;
clc;
set(groot,'defaultLineLineWidth',2.0)

T0 = 0;
Tpi = 1;

N = [4, 8, 16, 32, 64, 128, 256, 512];

% vector to store the errors
errors = zeros(1,length(N));

% to display exact solution
M = N(length(N)); % use detail as fine as largest input
X = zeros(M,1);
X(1) = T0;
for i = 2:M
    X(i) = i*pi/M;
end


subplot(2,1,1);
plot(X,f_heat_eq_exact(X))
hold on

for i = 1:length(N)
    [z, e] = f_heat_eq_est(N(i));
    % store error
    errors(i) = e;
%{
    % add graph to plot
    Xt = zeros(1, N(i)+1);
    Xt(1) = 0;
    for k = 2:N(i)
        Xt(k) = k*pi/N(i);
    end
    zt = [T0 transpose(z) Tpi];
    plot(Xt, zt)
%}
    Q = length(z);
    Xt = zeros(Q, 1);

    for k = 1:Q
        Xt(k) = k*pi/N(i);
    end
    plot([0;Xt;pi], [T0;z;Tpi])
    
end

legend("exact solution", "N = 4", "N = 8", "N = 16", "N = 32", ...
    "N = 64", "N = 128", "N = 264", "N = 512");
legend('Location','best')
xlabel("Position on heated rod")
ylabel("Temperature")
title("Numerical approximations of a differential heat equation")

subplot(2,1,2);
loglog(N, errors)
xlabel("Number of iterations")
ylabel("Error")
title("Error of approximations")
hold off


