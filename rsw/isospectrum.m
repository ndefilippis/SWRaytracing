function [FK] = isospectrum(Fkl);

% [FK] = ISOSPECTRUM(Fkl) 
%     Azimuthally integrate spectral field Fkl:  
%     FK = \int_0^{2\pi} Fkl K d\theta
%
%     Fkl assumed to be the result of fftshift(fft2(F)), with size(F)
%     = [nx,nx], nx a power of 2, followed by F = F(2:end,2:end), to
%     remove dummy values.  So size(Fkl) = [nx-1,nx-1].

[nkx,nky] = size(Fkl);
kmax = (nkx-1)/2;

[kx_,ky_] = ndgrid(-kmax:kmax,-kmax:kmax);
K_ = sqrt(kx_.^2+ky_.^2);
FK = zeros([kmax 1]);

for K = 1:kmax
    mask = (floor(K_+.5)-K)==0;
    FK(K) = sum(mask(:).*Fkl(:));
end 

