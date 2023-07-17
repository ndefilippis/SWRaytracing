
f = 3;   
C = 2;
time = linspace(0,10,100);

k0 = 20;    % x-wavenumber of carrier wave
a = .1;     % Amplitude of wavepacket
%sx = 10/L;  % Inverse widths of wavepacket in x and y
%sy = 10/L;
sx = 0;
sy = 0;

nx = 256;   % Make it a power of 2
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx) - L/2;
dx = 2*pi/nx;     % L_domain = 2*pi
[x_,y_] = ndgrid(x,x);

w0 = sqrt(f^2+C^2*k0^2);
cp = w0/k0;
cg = C/cp;

x0 = 0;
y0 = 0;
henv = a*exp(-((x_-x0)*sx).^2-((y_-y0)*sy).^2);

clear Ui
Ui(:,:,1) = w0/k0*henv.*cos(k0*x_);
Ui(:,:,2) = f/k0*henv.*sin(k0*x_);
Ui(:,:,3) = henv.*cos(k0*x_);

[U] = lsw(Ui,f,C,time);

clear M
figure
pcolor(Ui(:,:,3)'), shading interp, axis image, colorbar
ca = caxis;
for it = 1:length(time)
    pcolor(U(:,:,3,it)'), shading interp, axis image, caxis(ca), colorbar
    M(it) = getframe;
end
%mplay(M)
%gifmovie(M,'foo',.05)


 
% PV q = (v_x + f)/(1+h)
%kf = [0:kmax -kmax-1:-1]';
%q = (real(ifft(i*repmat(kf,[1 size(U,3)]).*fft(squeeze(U(:,2,:)))))+f)./(1+squeeze(U(:,3,:)));