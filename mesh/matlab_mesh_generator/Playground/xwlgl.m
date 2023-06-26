function [x,w] = xwlgl(np,a,b) 
%XWLGL Computes nodes and weights of the Legendre-Gauss-Lobatto quadrature formula.
%
%    [x,w]=xwlgl(np) returns the np weigths and nodes 
%    of the corresponding Legendre Gauss-Lobatto quadrature 
%    formula in the reference interval [-1,1].
%
%    [x,w]=xwlgl(np,a,b) returns the np weigths and the nodes 
%    of the corresponding Legendre Gauss-Lobatto quadrature 
%    formula in the  interval [a,b].
%
% Input: np = number of nodes
%        a, b = extrema of the interval
%
% Output: x(np,1) = LGL nodes  (CHQZ2, (2.3.12), pag. 76)
%                   the set is given by (CHQZ2, pag. 92) as well:
%               {-1} U {zeros of Jacobi polynomial P_{N-1}^{(1,1)}} U {1}
%         w(np,1) = LGL weigths (CHQZ2, (2.3.12), pag. 76)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


if np <= 1 
  es='The number of the quadrature nodes should be greater than 1';
  disp(es); x = []; w = [];
  return 
  
elseif np==2
    x=[-1;1];w=[1;1];
    
else
x=zeros(np,1);
w=zeros(np,1);
n=np-1;
x(1)=-1;x(np)=1;
x1=jacobi_roots(n-1,1,1);
x(2:n)=x1;
coef=2/(n*np);
w=coef./(pnleg(x,n).^2);
end

if nargin ==3 
   bma=(b-a)*.5;
   bpa=(b+a)*.5;
   x=bma*x+bpa;
   w=w*bma;
end

return
