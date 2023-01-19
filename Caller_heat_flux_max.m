clear;
clc;
format long;

% REMEMBER TO TAKE MINUS THE FUNCTION!!

tol = 1e-6;
iter_max = 1000;
needed_N = NaN;
final_p = NaN;
q_max = NaN;
a = 0; % ????
b = 2; % ????
% take first iteration for comparison, qc = q_current, pc = P_current
[pc, qc] = f_min_golden(@(p)f_heat_flux(p,3), a, b, tol);


for i = 4:iter_max

    [p, q] = f_min_golden(@(p)f_heat_flux(p,i), a, b, tol);
    rel_err = abs((p-pc) / pc);
    if (rel_err < tol)
        needed_N = i;
        final_p = p;
        q_max = - q;
        return
    end
    % if not, update your pc and qc
    pc = p;
    qc = q;
end

% if you make it this far, something's wrong
disp("Maximum not found. Consider adjusting initial bracket for p.")

