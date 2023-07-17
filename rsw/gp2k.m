function prodk = gp2k(prodg,icalphak)
    
    % Transform grid field, holding product of fields, 
    % with imaginary part holding field
    % shifted by dx/2, to k-space, and average result
        

    nx=size(prodg,1);
    kmax=nx/2-1;
    
    Wk  = fftshift(fft2(prodg))/nx^2;
    
    % Extract spectral products on grid and shifted grid, and average.
    
    Wk_up = Wk(2:end,kmax+2:end);
    Wk_dn = rot90(rot90((conj(Wk(2:end,2:kmax+2)))));
    prodk = ((1 - icalphak).*Wk_up + (1 + icalphak).*Wk_dn)/4;

    return

