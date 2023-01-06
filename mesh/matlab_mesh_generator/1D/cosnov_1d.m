function [nov]=cosnov_1d(npdx,ne);
% COSNOV_1D : Constructs the 1D (local mesh ---> global mesh) map
% 
% [nov]=cosnov_1d(npdx,ne,nov) computes the 2-indeces array nov:
%       nov(i,ie) = is the index, with respect to global numbering,
%                   associated to node i of element ie.
%
% Input: npdx = number of nodes in one element (the same in every element)
%        ne   = number of elements
%        nov  = zeros(max(npdx),ne)
% Output: nov = 2-indeces array of local to global map.
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$
nov = zeros(npdx, ne);
nov(1:npdx,1)=(1:npdx)';
k=npdx(1);
for ie=2:ne
k1=npdx-1;
nov(1:npdx,ie)=(k:k+k1)'; k=k+k1;
end
