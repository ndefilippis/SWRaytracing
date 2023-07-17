function fkf = fulspecs(fk);

    %     Assumes fk contains upper-half plane of spectral field, 
    %     and specifies lower half plane by conjugate 
    %     symmetry.  Input fk should be (1:2*kmax+1,1:kmax+1), kmax = 2^n-1.
    %     fkf is padded with zeros on (1,:) and (:,1), as expected
    %     by fftshift.  Grid resolution will be 2^(n+1) x 2^(n+1).  

    
    [nx,nky] = size(fk);
    kmax = nky-1;

    fkf = zeros(nx,nx);
    fkf(1:end,1:kmax+1) = fk;
    fkf(kmax+3:end,1) = conj(fkf(kmax+1:-1:2,1));
    fkf(1,kmax+3:end) = conj(fkf(1,kmax+1:-1:2));
    fkf(2:kmax+1,kmax+3:end) = conj(fkf(end:-1:kmax+3,kmax+1:-1:2));
    fkf(kmax+3:end,kmax+3:end) = conj(fkf(kmax+1:-1:2,kmax+1:-1:2)); 
        
    return
    
