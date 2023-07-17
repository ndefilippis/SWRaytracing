function [KEw,PEw,KEg,PEg,kes,pes,kesw,pesw,kesg,pesg,Um,etawk,ugk,vgk] = wavevortdecomp(U,f,c)

%  [....] = wavevortdecomp(U,f,c)  Decompose fields into two wave
%  parts and a geostrophic part, using linear theory.

[nx,ny,nf,nt] = size(U);  % nf = 3
kmax = nx/2-1;
[kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);

K2_   = kx_.^2 + ky_.^2;
sig2_ = f^2 + c^2*K2_;
iK2_ = 1./K2_;
iK2_(kmax+1,1) = 0;  % remove NaN

clear kes pes kesw pesw KEw PEw KEg PEg etawk
for j=1:nt
    uk   = g2k(U(:,:,1,j));
    vk   = g2k(U(:,:,2,j));
    etak = g2k(U(:,:,3,j));
    
    Um(j,1)=uk(kmax+1,1);
    Um(j,2)=vk(kmax+1,1);

    kes(:,j) = iso_spectra(uk.*conj(uk)+vk.*conj(vk));
    pes(:,j) = c^2*iso_spectra(abs(etak).^2);
    
    deltak = 1i*(kx_.*uk + ky_.*vk);
    zetak  = 1i*(kx_.*vk - ky_.*uk);
    
    etawk(:,:,j) = (f*zetak + c^2*K2_.*etak)./sig2_;
    %zetaw = f*etaw;
    
    KEwk = (abs(deltak).^2+f^2*abs(etawk(:,:,j)).^2).*iK2_;
    KEw(j) = sum(sum(KEwk));
    PEwk = c^2*abs(etawk(:,:,j)).^2;
    PEw(j) = sum(sum(PEwk));

    kesw(:,j) = iso_spectra(KEwk);
    pesw(:,j) = iso_spectra(PEwk);
    
    etagk = (f*etak-zetak).*f./sig2_;
    zetagk = -c^2/f*etagk.*K2_;
    
    ugk = -1i*ky_.*(c^2/f*etagk);
    vgk = 1i*kx_.*(c^2/f*etagk);  
    
    KEgk = abs(zetagk).^2.*iK2_;
    PEgk = c^2*abs(etagk).^2;
    KEg(j) = sum(KEgk(:));
    PEg(j) = sum(PEgk(:));

    kesg(:,j) = iso_spectra(KEgk);
    pesg(:,j) = iso_spectra(PEgk);

end



%fg  = nx^2*real(ifft2(ifftshift(fk)));
