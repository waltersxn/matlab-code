function z = f_optic_t(X)
% evaluates the time of optic ray transition defined in 
% problem 1 of HW5


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

% initialize your output
z = 0;

% X is the vector of interior points,
% we evaluate from the first to the penultimate
for i = 1:(length(X) + 1)
    
    z = z + sqrt( (y(i+1) - y(i))^2 + (x(i+1) - x(i))^2 )  / v(i);

end

end