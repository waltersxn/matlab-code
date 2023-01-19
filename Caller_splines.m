clear;
clc;


for n = 2:4

X = zeros(n+1,1);

h = 2*pi/n;

for j = 1:n+1
    X(j) = (j-1)*h;
end

Y = cos(X);

% return complete, not-a_knot, your spline for this n

x = linspace(0,2*pi,101);
% subplot
C = f_Qspliner(X, Y, 0);

y = f_splinerator(C, X, x);

spline()

plot(x,y, x, cos(x))
hold on

end
hold off