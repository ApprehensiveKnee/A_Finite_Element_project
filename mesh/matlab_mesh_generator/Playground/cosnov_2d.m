function [nov]=cosnov_2d(npdx,nex,npdy,ney);
% COSNOV_2D : Constructs the 2D (local mesh ---> global mesh) map
%
% [nov]=cosnov_2d(npdx,nex,npdy,ney) computes the 2-index array nov:
%       nov(i,ie) = is the index, with respect to global ordering,
%                   associated to node i of element ie.
%
%       __________________________
%       |      |      |     |     |
%       |  3   |  6   |  9  | 12  |      Omega and spectral elements 
%       |      |      |     |     |      ordering
%       __________________________
%       |      |      |     |     |
%       |  2   |  5   |  8  | 11  |    
%       |      |      |     |     |
%       __________________________
%       |      |      |     |     |
%       |  1   |  4   |  7  | 10  |
%       |      |      |     |     |
%       __________________________
%
%
% Input: npdx = number of nodes along x-direction inside 
%               one element (the same in every element)
%        nex  = number of elements along x-direction
%        npdy = number of nodes along y-direction inside
%               one element (the same in every element)
%        ney  = number of elements along y-direction
%
% Output: nov = 2-index array of local to global map, 
%               size(nov)=[max(npdx*npdy),ne]
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

ldnov=npdx*npdy;
ne=nex*ney;
nov=zeros(ldnov,ne);

% element 1

nov(:,1)=(1:ldnov)';

% elements first column

k=ldnov-npdx;
for ie=2:ney
nov(:,ie)=(k+1:k+ldnov)';
k=k+ldnov-npdx;
end
kmax=k+npdx;

% other columns

nxm1=npdx-1;
for iex=2:nex

%  other rows, bottom elements

ie=(iex-1)*ney+1;
for j=1:npdy
k=(j-1)*npdx;
nov(k+1,ie)=nov(j*npdx,ie-ney);
nov(k+2:k+npdx,ie)=(kmax+1:kmax+nxm1)';
kmax=kmax+nxm1;
end

% other elements

for iey=2:ney
ie=(iex-1)*ney+iey;

% first row
nov(1:npdx,ie)=nov(ldnov-npdx+1:ldnov,ie-1);
for j=2:npdy
k=(j-1)*npdx;
nov(k+1,ie)=nov(j*npdx,ie-ney);
nov(k+2:k+npdx,ie)=(kmax+1:kmax+nxm1)';
kmax=kmax+nxm1;
end

end
end
