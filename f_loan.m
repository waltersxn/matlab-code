function q = f_loan(rate)
  % LOAN evaluates normalized balance for interest rate r, amount a,
  % payment p, after n years
  amount = 100000;
  payment = 10000;
  n = 20;


  q = (1 + rate).^n - (payment/amount)*((1+rate).^n-1)/rate;

  loan(rate == 0) = amount - (n*payment);
end
  % as discussed in class, this formula will have divide by zero when r=0!
  % so you will need to modify it to handle that special case. keep in mind
  % that this function is vectorized, so using a test like
  % if r==0, f=...; end
  % will not work. instead you can use the syntax
  % f(r==0) = ...;
  % to overwrite only the entry in f for which r=0
  % (or you can loop over all entries and perform the test for each entry
