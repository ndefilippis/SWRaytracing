function fg = k2g(fk)

% fg = k2g(fk)

nx = size(fk,1)+1;
fg  = nx^2*ifft2(ifftshift(fulspec(fk)));