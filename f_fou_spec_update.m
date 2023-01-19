function u_t = f_fou_spec_update(t, N, delta_t, u_0)
% t = desired time of result
% delta_t = size of time step
% N = number of grid points

% check:: u_0 is a vector of size N
if ( length(u_0) ~= N)
    disp("Your input vector length must match your number of grid points.")
    return
end

% number of time steps
if (delta_t > 0 && t > 0)
    time_steps = ceil(t / delta_t);
else
    disp("Time step, delta_t must be positive");
    return
end

% get fourier spectral differentiation matrix
D = f_four_spec_diff_mat(N);
% and identity matrix
I = eye(N);

% find the stepping matrix based on our calculation
% u(t_n) = (delta_t * D + I_N)u(t_{n-1})
A = D.*delta_t + I;
[L, U] = lu(A);

% solve with iterations of back/forward substitution
u_t = u_0;

for i = 1:time_steps
    y = L \ u_t;
    u_t = U \ y;
end


end