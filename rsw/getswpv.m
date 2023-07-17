function [q,zeta,qlin] = getswpv(U,f);

%  [q,zeta,lin] = getswpv(U,f)  U is state vector from swk output. 

[nx,ny,nf,nt] = size(U);  % nf = 3
kmax = nx/2 - 1;  % 2^n - 1
nkx = 2*kmax+1;   % -kmax:kmax
nky = kmax+1;     % 0:kmax
[ikx_,iky_] = ndgrid(-kmax:kmax,0:kmax);
K_   = sqrt(ikx_.^2 + iky_.^2);
ikx_ = sqrt(-1)*ikx_; 
iky_ = sqrt(-1)*iky_;

q = zeros(nx,ny,nt);
zeta = q;
for j = 1:nt
    zeta(:,:,j) = k2g(-iky_.*g2k(U(:,:,1,j))+ikx_.*g2k(U(:,:,2,j)));
    q(:,:,j) = (zeta(:,:,j)+f)./(1+U(:,:,3,j));
    qlin(:,:,j) = zeta(:,:,j)-f*U(:,:,3,j);
end

 