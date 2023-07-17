function [dfdx] = dxk(f);

% Compute df/dx using spectral transform

f = squeeze(f);
i = sqrt(-1);
kmax = size(f,1)/2-1;
kf = repmat([0:kmax -kmax-1:-1]',[1 size(f,2)]);
dfdx = real(ifft(i*kf.*fft(f)));
