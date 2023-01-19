clear;
clc;
set(groot,'defaultLineLineWidth',2.0)


% these are all constants and do not change
tol = 1e-6;
pw = 10^5;
rhow = 10^3;
Req = 1;
tmax = 1;
R0 = 0.1;
dR0 = 0;
gamma = 1.4;
Y0 = [R0, dR0];

% time step and count for non-adaptive
dt = 1/16000;
N = tmax/dt;


% crude (non-adaptive) implementation
Ycurr = Y0;
Ynon = zeros(N+1,1);
Ynon(1) = Ycurr(1);
for n = 1:N

    K1 = -f_Dy(Ycurr);
    K2 = -f_Dy(Ycurr +0.5*dt*K1);
    K3 = -f_Dy(Ycurr +0.5*dt*K2);
    K4 = -f_Dy(Ycurr +dt*K3);

    Ycurr = Ycurr + dt/6 * (K1 + 2*K2 + 2*K3 + K4);
    Ynon(n+1) = Ycurr(1);
end

Xnon = linspace(0, tmax, length(Ynon));

% adaptive timestep implementation
Dt_ad = 1e-7; % the initial time step value
t = 0; % the current time
T = t; % a vector for times, also use for step count
% steps = 0; % counter for #steps
Ycur_ad = Y0;
% solution vector for adaptive method 
Yad = Y0(1); % will become the solution vector
% vector of local error estimates
lee = 0;
% vector of Dt's 
Tdt = Dt_ad;
% remember to plot these starting at index 2

while (t < tmax)
    
    K1ad = f_Dy(Ycur_ad);
    K2ad = f_Dy(Ycur_ad + 0.5*Dt_ad*K1ad);
    K3ad = f_Dy(Ycur_ad - Dt_ad*K1ad + 2*Dt_ad*K2ad);
    Ysimp = Ycur_ad + 0.5*Dt_ad*(K1ad + K3ad);
    Yhat = Ycur_ad + Dt_ad*(K1ad + 4*K2ad + K3ad) / 6;
    l2dt = norm(Ysimp - Yhat);

    if (l2dt > tol)
        Dt_ad = 0.9*Dt_ad*(tol/l2dt)^(1/3);
    else
        Ycur_ad = Yhat;
        lee = [lee, l2dt];
        Yad = [Yad, Ycur_ad(1)];
        Tdt = [Tdt, Dt_ad];
        t = t + Dt_ad;
        T = [T,t];
        Dt_ad = 0.9*Dt_ad*(tol/l2dt)^(1/3);
    end
end
Tdt = Tdt(2:length(Tdt));
lee = lee(2:length(lee));

subplot(2,2,1)
plot(Xnon,Ynon, T, Yad);
hold on
title("Solution plotting of u(x)")
legend("Non-adaptive, N = 16000", "Adaptive")
xlabel("x")
ylabel("u(x)")

subplot(2,2,2)
plot(lee);
title("Local error estimate")
subplot(2,2,3)
plot(Tdt)
title("Size of time step vs. time step #")
xlabel("# of time steps")
ylabel("Size of dt")

% the function used to calculate matrix-free differentiation
% does not depend on dt so is same for both methods
function Dy = f_Dy(Y)
% repeat all the same constants
pw = 10^5;
rhow = 10^3;
Req = 1;
tmax = 1;
R0 = 0.1;
dR0 = 0;
gamma = 1.4;

Dy = zeros(2,1);
Dy(1) = Y(2);
Dy(2) = ((pw/rhow)*(Req/Y(1))^(3*gamma) - (pw/rhow) ...
    - (3/2)*(Y(2))^2) / Y(1);

end


