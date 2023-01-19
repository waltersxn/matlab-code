% CHEB  compute D = differentiation matrix, x = Chebyshev grid
% from https://people.maths.ox.ac.uk/trefethen/spectral.html
  function [D,x] = f_cheb(N)
  if N==0, D=0; x=1; return, end
  x = cos(pi*(0:N)/N)'; 
  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
  X = repmat(x,1,N+1);
  dX = X-X';                  
  D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
  D  = D - diag(sum(D'));                 % diagonal entries
  % above is standard, with x(1)=1 and x(N+1)=-1
  % (points are reversed compared to what might be expected)
  % below I'm changing to the more familiar order,
  % though you would not do this in an efficient implementation
  x = flipud(x);
  D = flipud(fliplr(D));