function q = f_heat_flux(p, n)
% given an input parameter p and a number of partitions of the heat rod
% assumed to have lenght pi, delivers the heat flux q at the 0 end
% of the rod. This is provided as an estimate of the first derivative
% of T at point 0, where T is the heat equation on the rod.


    % first, still need to get the vector of Ti's
    T0 = 0;
    Tn = -exp(p);
    h = pi/n;

    % make your vector of x_i's (partition the rod)
    X = zeros(n-1, 1);
    for i = 1:n-1
        X(i) = i*pi/n;
    end
    
    % populate the matrix
    % this time, putting h^2 in the b vector

    A = zeros(n-1, n-1);
    % edge cases
    A(1,1) = -2;
    A(1,2) = 1;
    A(n-1, n-2) = 1;
    A(n-1,n-1) = -2;
    % diagonal and sub/superdiagonals
    for i = 2:n-2
        A(i,i) = -2;
        A(i,i-1) = 1;
        A(i,i+1) = 1;
    end

    % populate b vector NOTE:: using -psin(x)
    b = (h^2)*(-p)*sin(X);
    b(1) = b(1) + T0;
    b(n-1) = b(n-1) - Tn;
    
    % solve for the vector of Ti's
    
    Z = A\b;

    % use the first two elts to approximate T'(0)
    % recall that h is the "step" size
    % first order finite difference
    q = - (Z(2) - Z(1)) / h; % setting minus to return the minimizing func

end