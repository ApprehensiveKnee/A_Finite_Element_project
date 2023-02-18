
function [A]=stiff_2d_sp(wx,dx,jacx,wy,dy,jacy);
% STIFF_2D_SP Computes 2D local SEM matrix associated to (nabla(phi_j), nabla(phi_i))_N
%
%     (nabla(phi_j), nabla(phi_i))_N 
%
%
%    [A]=stiff_2d_sp(wx,dx,jacx,wy,dy,jacy);
%        produces the matrix
%        A_{ij}=(nabla(phi_j), nabla(phi_i))_N 
%        of size
%        (mn,mn) where mn=npdx*npdy is the local number of d.o.f.
%
% Input : 
%         wx = npdx LGL weigths in [-1,1],
%            (produced by calling [x,wx]=xwlgl(npdx))
%         dx =first derivative LGL matrix (by calling dx=derlgl(x,npdx))
%         jacx = array (length(jacx)=ne); jacx(ie)= (x_V2_ie-x_V1_ie)/2
%         wy = npdy LGL weigths in [-1,1],
%            (produced by calling [y,wy]=xwlgl(npdy))
%         dy =first derivative LGL matrix (by calling dy=derlgl(y,npdy))
%         jacy = array (length(jacy)=ne); jacy(ie)= (y_V3_ie-y_V1_ie)/2
%         ww = column array with local weigths, length: npdx*npdy
%
% Output: A = matrix (npdx*npdy,npdx*npdy)
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

npdx=length(wx);
npdy=length(wy);
jacyx=jacy/jacx;
jacxy=jacx/jacy;
mn=npdx*npdy;
[wx1,wy1]=meshgrid(wx,wy); ww=wx1.*wy1; ww=ww'; ww=ww(:); clear wx1 wy1;

A=sparse(mn,mn);A=0;
B=sparse(mn,mn);B=0;
for ki=1:npdy
    inde=((ki-1)*npdx+1:ki*npdx);
    A(inde,inde)=dx'*(diag(ww(inde))*dx)*jacyx;
end

for i=1:npdx
    inde=(i:npdx:mn);
    B(inde,inde)=dy'*(diag(ww(inde))*dy)*jacxy;
end
A=A+B;

