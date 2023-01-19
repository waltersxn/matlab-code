function Jt = f_gradPhi_jac(X)
% returns a Jacobian for f(x) = gradient(phi(x)) from HW4 prob 1e)
% Assumed: vector X is length 2

if (length(X) ~= 2)
    disp("Phi_jac:: input X must be lenght 2!!");
    return
end

% first, your vector of ti's and bi's

t = [0, 1, 2];
b = [exp(0.1); exp(0.9); exp(2)];

% build and populate the residual vector for phi(x) = 1/2(||r(x)||^2)
r = zeros(length(t), 1);
for i = 1:length(t)
    r(i) = X(1)*exp(t(i)*X(2)) - b(i);
end

% build and populate the gradient of r

del_r = zeros(length(t), 2);
for i = 1:length(t)
    del_r(i,1) = exp(t(i)*X(2));
    del_r(i,2) = t(i)*X(1)*exp(t(i)*X(2));
end

% this should equal phi_grad:: verify
%Phi = transpose(del_r)*r % verified!!

% this is the matrix that Gauss-Newton uses
L = transpose(del_r) * del_r;

% this will be the "rest" of the full Hessian
d2_r1 = zeros(2, 2);

d2_r2 = zeros(2, 2);
d2_r2(1,2) = exp(X(2));
d2_r2(2,1) = exp(X(2));
d2_r2(2,2) = X(1)*exp(X(2));

d2_r3 = zeros(2, 2);
d2_r3(1,2) = 2*exp(2*X(2));
d2_r3(2,1) = d2_r3(1,2);
d2_r3(2,2) = 4*X(1)*exp(2*X(2));

M = r(1)*d2_r1 + r(2)*d2_r2 + r(3)*d2_r3;

Jt = L + M;

% old version of J that is no good
%for i = 1:3
 %   J(1,1) = J(1,1) + exp(2*X(2)*t(i));
  %  J(1,2) = J(1,2) + 2*t(i)*X(1)*exp(2*X(2)*t(i)) - b(i)*t(i)*exp(X(2)*t(i));
   % J(2,1) = J(2,1) + 2*t(i)*X(1)*exp(2*X(2)*t(i)) - b(i)*t(i)*exp(X(2)*t(i));
    %J(2,2) = J(2,2) + 2*(t(i))^2*(X(2))^2*exp(2*X(2)*t(i)) ...
     %   - b(i)*(t(i))^2*X(1)*exp(X(2)*t(i));
%end



end