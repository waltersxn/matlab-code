
%r = 0.01:0.01:1;
%rr = f(r);
%rdr = der_f(r);
%plot(r,rr, 'b', r, rdr, 'r')
err = 1e-6;
start_at = 0.5;

[zer, Points] = newton(@loan, @loan_der, start_at, err);


t = err:length(Points);
normt = t./(length(Points)-1);


semilogy(normt, Points, 'g')
% hold on
% semilogy(your other stuff)



    

