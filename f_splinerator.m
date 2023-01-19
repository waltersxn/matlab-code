function y = f_splinerator(C, X, x)
% evaluates x with the interpolated function constructed from
% input C: the nx3 matrix of coefficients of quadratic eqns
% X: the vector of points for which C was evaluated
% x: the point at which the function should be evaluated

y = zeros(length(x), 1);

for k = 1:length(x) 
    z = x(k); % for each value in your input vector
    for i = 1:(length(X)-1) 
        if (X(i) <= z) && (z <= X(i+1)) % find the interval cont. z, evaluate
            y(k) = C(i,1) + C(i,2).*(z-X(i)) + C(i,3).*(z-X(i)).^2;
            break
        end
    end

end

end