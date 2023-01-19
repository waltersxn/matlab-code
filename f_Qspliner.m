function C = f_Qspliner(X, Y, yp0)
% inputing X: the vector of xi's; Y: the vector of yi's
% yp0: the first derivative at x0

% output, an nx3 matrix of coefficients for the piecewise 
% quadratic interpolation of the function that assigned the xi's to yi's

if (length(X) ~= length(Y))
    disp("Vectors X and Y must have equal length")
    return
end

if (length(X) < 3)
    disp("Must have at least 3 points")
    return
end


h = X(2) - X(1); % this is fixed for all xi
n = length(X) - 1;
yp = yp0;

C = zeros(n,3); % the matrix to load our coefficients
U = zeros(3,3); % the matrix to solve updates with

U(1,1) = 1;
U(2,1) = 1;
U(2,2) = h;
U(2,3) = h^2;
U(3,2) = 1;

for i = 1:n
    
    b = [Y(i); Y(i+1); yp]; % set/update your output vector
    x = U \ b; % get these three coefficients
    % update the Coefficient matrix
    C(i,1) = x(1);
    C(i,2) = x(2);
    C(i,3) = x(3);
    
    yp = x(2) + 2*h*x(3); % update the derivative approximation

end

end