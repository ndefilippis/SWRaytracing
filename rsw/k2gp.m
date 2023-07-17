function fg = k2gp(fk,damask,ialphakf)

   % Transform to grid space, packing fields shifted by dx/2 into
   % imaginary part
    
   nx = size(fk,1)+1;
   fkt = fulspec(damask.*fk).*(1 + ialphakf);
   fg  = nx^2*ifft2(ifftshift(fkt));
