clear;
clc;

a1 = pi/3;
a2 = 7*pi/3;
n_vec = [2; 3; 5; 10; 20];
N = length(n_vec);

T_1out = transpose(zeros(1, N));
T_1error = transpose(zeros(1, N));

T_2out = transpose(zeros(1, N));
T_2error = transpose(zeros(1, N));

for i = 1:N
    [T_1out(i), T_1error(i)] = f_taylor_sin(a1, n_vec);
end

for i = 1:N
    [T_2out(i), T_2error(i)] = f_taylor_sin(a2, n_vec);
end

%[T_est, T_error] = f_taylor_sin(a1, n);


%results = table(n_vec,T_est,T_error);
