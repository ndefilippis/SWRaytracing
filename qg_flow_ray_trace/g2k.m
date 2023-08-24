function fk = g2k(fg)
    
% fk = g2k(fg)
    
nx = size(fg,1);
kmax = nx/2 - 1;
    
fkt = fftshift(fft2(fg))/nx^2;
fk = fkt(2:end,kmax+2:end);      

