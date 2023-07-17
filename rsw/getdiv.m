function delta = getdiv(U,f);

%  delta = getdiv(U)  U is state vector from swk output. 

[nx,ny,nf,nt] = size(U);  % nf = 3
kmax = nx/2 - 1;  % 2^n - 1
nkx = 2*kmax+1;   % -kmax:kmax
nky = kmax+1;     % 0:kmax
[ikx_,iky_] = ndgrid(-kmax:kmax,0:kmax);
%K_   = sqrt(ikx_.^2 + iky_.^2);
ikx_ = sqrt(-1)*ikx_; 
iky_ = sqrt(-1)*iky_;

delta = zeros(nx,ny,nt);
for j = 1:nt
    delta(:,:,j) = k2g(ikx_.*g2k(U(:,:,1,j))+iky_.*g2k(U(:,:,2,j)));
end

