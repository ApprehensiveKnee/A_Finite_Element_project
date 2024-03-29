function [err_l2]=normal2_2d(fdq,nq,errtype,x,wx,xx,jacx,y,wy,yy,jacy,...
xy,ww,nov,un,u,uex);
% NORMAL2_2D   Computes L2-norm in 2D
%
%  [err_l2]=normal2_2d(fdq,nq,errtype,x,wx,xx,jacx,y,wy,yy,jacy,...
%                     xy,ww,nov,un,u);
%
% Input : fdq = 0  uses Legendre Gauss quadrature formulas with nq nodes
%               in each element (exactness degree = 2*nq+1)
%             = 1  uses Legendre Gauss Lobatto quadrature formulas with npdx/y
%                  quadrature nodes in each element.
%                Quadrature nodes are the nodes of the mesh.
%         nq = nodes (in each element and along each direction) 
%              for GL quadrature formulas. Not used if  fdq == 1
%         errtype = 0 for absolute error ||u-u_ex||
%                   1 for relative error ||u-u_ex||/||u_ex||
%         x = column array  with npdx LGL nodes in [-1,1]
%         wx= column array  with npdx LGL weigths in [-1,1]
%         dx= derivative matrix produced with derlgl
%         xx = 2-indexes array of size (4,ne) 
%            xx(1:4,ie)=[x_V1_ie;x_V2_ie;x_V3_ie;x_V4_ie]
%            (ne=nex*ney is the global number of spectral elements)
%         jacx = array (length(jacx)=ne); jacx(ie)= (x_V2_ie-x_V1_ie)/2
%         y = npdy LGL nodes in [-1,1], previously generated by xwlgl
%         wy = npdy LGL weigths in [-1,1], previously generated by xwlgl
%         dy= derivative matrix produced with derlgl
%         yy = 2-indexes array of size (4,ne):
%            yy(1:4,ie)=[y_V1_ie;y_V2_ie;y_V3_ie;y_V4_ie]
%         jacy = array (length(jacy)=ne); jacy(ie)= (y_V3_ie-y_V1_ie)/2
%         xy = column array with global mesh, length: noe=nov(npdx*npdy,ne)
%         ww = column array with global weigths, length: noe=nov(npdx*npdy,ne)
%              diag(ww) is the mass matrix associated to SEM discretization
%         nov = local to global map. 2-indexes array, size(nov)=[nov,ne]
%         un = column array with the numerical solution
%         u  = column array with the evaluation of exact solution
%         uex  = exact solution (uex=@(x,y)[uex(x,y)], with .*, .^, ./)
%
% Output: err_l2 = error in L2-norm
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$



npdx=length(x); npdy=length(y);
[ldnov,ne]=size(nov); noe=nov(ldnov,ne);
num=0; den=0;
err_l2=0;

if fdq==0

% Legendre Gauss nodes and weigths

[xg,wxg] = xwlg(nq,-1,1);
[yg,wyg] = xwlg(nq,-1,1);
[wxg1,wyg1]=meshgrid(wxg,wyg); wxyg=wxg1.*wyg1;wxyg=wxyg';wxyg=wxyg(:); 
clear wxg1 wyg1;

% Evaluation of Lagrange basis polynomials at quadrature nodes xg
%
[phix]= intlag_lgl (x,xg);
[phiy]= intlag_lgl (y,yg);

% Loop on spectral elements

for ie=1:ne
un_loc=un(nov(1:ldnov,ie));
[norma_loc1,norma_loc2]=normal2_ex_loc(errtype,nq,uex,un_loc,wxyg,...
 jacx(ie),xx(1:4,ie),xg,phix,...
 jacy(ie),yy(1:4,ie),yg,phiy);
num=num+norma_loc1;
den=den+norma_loc2;
end

elseif fdq==1
    
num=sum((u-un).^2.*ww);
den=sum(u.^2.*ww);

end


if errtype==0
    err_l2=sqrt(num);
elseif errtype==1
    if abs(den)>1.d-14; err_l2=sqrt(num/den); end
end

return

function [norma_loc1,norma_loc2]=normal2_ex_loc(errtype,nq,uex,...
 un,wxyg,jacx,xx,xg,phix,...
    jacy,yy,yg,phiy);

% High degree Legendre Gaussian  formulas to compute H1-norm error

[nq,npdx]=size(phix);
[nq,npdy]=size(phiy);

% mapping quadrature nodes on element ie

xxg=xg*jacx+(xx(2)+xx(1))*.5;
yyg=yg*jacy+(yy(3)+yy(1))*.5;
[xg1,yg1]=meshgrid(xxg,yyg); xg1=xg1'; yg1=yg1'; xyg=[xg1(:),yg1(:)]; 
clear xg1; clear yg1;


% evaluation of exact solution at quadrature nodes.

[U]=uex(xyg(:,1),xyg(:,2));


% evaluate numerical solution  at quadrature nodes.

un=reshape(un,npdx,npdy);

un_i=phix*(phiy*un')';

un_i=un_i(:);


% compute the sum

norma_loc1=sum((U-un_i).^2.*wxyg)*jacx*jacy;

if errtype==0
    norma_loc2=0;
else
    norma_loc2=sum(U.^2.*wxyg)*jacx*jacy;
end

return


