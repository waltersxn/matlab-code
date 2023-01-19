clear;
clc;

T = [12.420601220444135, 23.93446955965803, 25.819341664737628, 12.0];

[h, tvec] = tide_height_time();

%%%%%%%% first section %%%%%%%%%%%%%%%%
p = length(T);
A = zeros(length(h), 2*p + 1);
% note: for each odd 1 <= i < 2p, x_i = a_i * cos(phi_i)
% and x_{i+1} = a_i * sin(phi_i)
% such that A is linear wrt 
% X = [a_1*cos(phi_1), a_1*sin(phi_1), ... ]

% populate the matrix with coefficients of the linearization
for i = 1:length(h)
    for j = 1:2:(2*p-1)
        A(i,j) = cos(2*pi*tvec(i) / T( (j+1) / 2 ));
        A(i,j+1) = -sin(2*pi*tvec(i) / T((j+1)/ 2));
        A(i, 2*p + 1) = 1;
    end % column loop
end % row loop

% build the transpose and new b vector
B = transpose(A)*A;
b = transpose(A)*h;

% solving with Gaussian elimination
Xg = f_gauss_mat(B, b);

% build a vector whose elements are the reversal of the linearization
Yg = zeros(length(Xg), 1);
for j = 1:2:(2*p-1)
        Yg(j) = atan2(Xg(j+1), Xg(j));
        Yg(j+1) = Xg(j) / cos(Yg(j));
end
Yg(2*p + 1) = Xg(2*p +1);


% solve using backslash
Xb = A \ h;
Yb = zeros(length(Xb), 1);
for j = 1:2:(2*p-1)
    Yb(j) = atan2(Xb(j+1), Xb(j));
    Yb(j+1) = sqrt((Xb(j))^2 + (Xb(j+1))^2);
end
Yb(2*p + 1) = Xb(2*p +1);

%%%%%%%%%%%%%% second section %%%%%%%%%%%%%%%%%%%%
% now, for progressively added solutions:
% first, chop the original A-matrix and h-vector down to vertical length
range = 720;
num_constituents = 4;
hc = h(1:range);
tq = tvec(1:range);
Ac = A(1:range,:);
Errors = zeros(range, num_constituents); % store residual vectors

% using r istead of p so as not to mix up the variables
for r = 1:num_constituents
    Aq = Ac(:,1:(2*r+1)); % horizontal chop

    % find the x-solution vector
    Xq = Aq \ hc;

    % reverse solve for the phase and amplitude values
    Yq = zeros(length(Xq), 1);
    for j = 1:2:(2*r-1)
        Yq(j) = atan2(Xq(j+1), Xq(j)); % phase_i
        Yq(j+1) = sqrt((Xq(j))^2 + (Xq(j+1))^2); % amplitude_i
    end
    Yq(length(Yq)) = Xq(length(Xq));

    tide_out = tide_height(Yq, tq); % using the function defined below
    Errors(:,r) = tide_out - hc; % populate one column of errors
    subplot(3,2,r)
    plot(tq, tide_out, tq, hc);
    legend("Predicted height", "True height")
    title("Predicting tide heights using " + r + " tidal constituents")
end

Err_norms = zeros(num_constituents, 1);
for i = 1:num_constituents
    Err_norms(i) = norm(Errors(:,i));
end
Err_x = linspace(1,num_constituents,num_constituents);
subplot(3,2,5)
plot(Err_x, Err_norms);
title("L2 norm of residual (error) vector vs. p")

function h = tide_height(X, t)
% the function for tide height based on a vector of amplitudes and phases
% the input vector should:
%       *be arranged with phase_1, ampl_1, phase_2, ampl_2, ... height_0
%       *have length 3, 5, 7, or 9

% vector of constituents
T = [12.420601220444135, 23.93446955965803, 25.819341664737628, 12.0];

pf = length(X);


h = X(length(X)); % initialize h with baseline value

for k = 1:2:(pf-2)
    h = h + X(k+1)*cos(2*pi.*t ./ T((k+1) / 2) + X(k));
end

end

