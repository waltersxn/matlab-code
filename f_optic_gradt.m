function Z = f_optic_gradt(X)
%f_optic_grad returns a length-3 vector output for the gradient
% of phi, defined in HW5-1

if (length(X) ~= 3)
    disp("Input must have length 3")
    return
end

% the vector of velocities
v = [8,6,2,1];
% vector of y's
y = [0,2,3,6,10];
% augmented x vector for counting
x = [0, transpose(X), 1];


Z = zeros(3,1);

% this is dt/dx_2
for i = 2:4
    
    Z(i-1) = ((x(i) - x(i+1)) / ...
        (v(i) * sqrt( (y(i+1)-y(i))^2 + (x(i+1) - x(i))^2 ) )) ...
        ...
        + ((x(i) - x(i-1)) / ...
        (v(i-1) * sqrt( (y(i) - y(i-1))^2 + (x(i) - x(i-1))^2  ) ));
end

end