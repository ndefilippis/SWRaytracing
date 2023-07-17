nx = 256; % must be a power of 2
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx)' - L/2;
kmax = nx/2 - 1;  % 2^n - 1
nkx = 2*kmax+1;   % -kmax:kmax
dx = 2*pi/nx;     % L_domain = 2*pi

f = 1;
Cg = 1;
a0 = .01;

km = 20;
K = [-km:-1 1:km];

Uin = zeros(nx,3);
for k=K
    wp = sqrt(f^2+Cg^2*k^2);
    ak = a0; %a0*f/wp;
    phi = rand*2*pi;  % random phase
    Uin(:,1) = Uin(:,1) + ak*wp/k*cos(k*x+phi);
    Uin(:,2) = Uin(:,2) + ak*f/k*sin(k*x+phi);
    Uin(:,3) = Uin(:,3) + ak*cos(k*x+phi);
end

np = 20;
Xpin = .2*randn(np,1);
[U,ke,pe,time,Xp,t,ue] = sw1(Uin,x,f,Cg,10000,100,Xpin);

nt = size(U,3);
figure
for j=1:nt
    plot(x,U(:,3,j))
    axis([-pi pi -2*max(U(:,3,1)) 2*max(U(:,3,1))])
    M(j) = getframe;
end

figure
plot(t,Xp)

% Compute full Stokes drift from stable spectrum of waves:
% us = sum_k a_k^2*w_k/(2*k);
% where a_k are coefficients for FFT of eta (maybe with factor of
% 2.. think it through)
% Full k-omega spectrum of u...  is it all still exact waves?
% time-avg of low-freq parts of stokes?  
% Work this all out for spectrum of waves in 2D SW.
% See if exact solution in 1D can be gotten in Fourier space, or at
% least understood that way, since single wave ends up as stable
% spectrum with a few wave modes. 
% Show shock as higher-order asymptotic solution...

figure
uk0 = fft(U(:,1,1))/nx;  % u
ukm = fft(U(:,1,floor(nt/2)))/nx;
ukf = fft(U(:,1,end))/nx;
loglog(0:kmax,abs(uk0(1:kmax+1)))
hold
loglog(0:kmax,abs(ukm(1:kmax+1)))
loglog(0:kmax,abs(ukf(1:kmax+1)))


% PV q = (v_x + f)/(1+h)

q = (dif(squeeze(U(:,2,:)),1,1)/dx+f)./(1+squeeze(U(:,3,:)));
