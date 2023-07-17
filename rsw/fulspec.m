function fkf = fulspec(fk);

    %     Assumes fk contains upper-half plane of spectral field, 
    %     and specifies lower half plane by conjugate 
    %     symmetry.  Input fk should be (1:2*kmax+1,1:kmax+1), kmax = 2^n-1.
    %     fkf is padded with zeros on (1,:) and (:,1), as expected
    %     by fftshift.  Grid resolution will be 2^(n+1) x 2^(n+1).  

    
    [nkx,nky] = size(fk);
    nx = nkx+1;
    kmax = nky-1;
    
    fkf = zeros(nx,nx);
    fup = fk;
    fup(kmax:-1:1,1) = conj(fup(kmax+2:nkx,1));
    fdn = conj(fup(nkx:-1:1,nky:-1:2));
    fkf(2:nx,nky+1:nx) = fup;
    fkf(2:nx,2:nky) = fdn;
        
