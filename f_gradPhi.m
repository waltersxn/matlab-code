function Y = f_gradPhi(X)
% returns the vector output of the gradient of 
% Assumed: vector X is length 2

if (length(X) ~= 2)
    disp("Phi_jac:: input X must be lenght 2!!");
    return
end

% first, your vector of ti's

t = [0, 1, 2];
b = [exp(0.1); exp(0.9); exp(2)];
Y = zeros(2,1);

% calculated the gradient_phi in a different way
for i = 1:3
    Y(1) = Y(1) + X(1)*exp(2*X(2)*t(i)) - b(i)*exp(X(2)*t(i));
    Y(2) = Y(2) + t(i)*(X(1))^2*exp(2*X(2)*t(i)) ...
        - b(i)*t(i)*X(1)*exp(X(2)*t(i));
end
%disp("This is grad_phi");
%Y

end