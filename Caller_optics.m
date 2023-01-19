




clear;
clc;

y = [0,2,3,6,10];



% your initial guess
xcurr = [0;0;0];

tol = 1e-6;
tol2 = 1e-6;

zsd = f_steepest_dec(@f_optic_t, @f_optic_gradt, xcurr, tol);

x = [0, transpose(zsd), 1];

zbfgs = f_BFGS(@f_optic_t, @f_optic_gradt, xcurr, tol2);

plot(x,y);








