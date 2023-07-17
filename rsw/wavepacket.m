
nx = 256; % must be a power of 2
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx) - L/2;
dx = 2*pi/nx;     % L_domain = 2*pi

kmax = nx/2 - 1;  % 2^n - 1
k_ = 0:kmax;

f = 1;
C = 1; % sqrt(gamma)
a = .01;  
k0 = 6;
w0 = sqrt(f^2+C^2*k0^2);
cp = w0/k0;
cg = C/cp;

x0 = pi/6;
henv = exp(-(10*(x-x0)/L).^2);

clear Uin
Uin(:,1)=a*w0/k0*henv.*cos(k0*x);
Uin(:,2)=a*f/k0*henv.*sin(k0*x);
Uin(:,3)=a*henv.*cos(k0*x);

np = 25;
Xpin = .3*randn(np,1);
[U,ke,pe,time,Xp,t,ue] = sw1(Uin,x,f,C,10000,100,Xpin);

mfig=figure('position',[700 600 700 400]);
clear M
for j=1:size(U,3)
    xs = time(j)*c1;
    js = floor(xs/dx);
    subplot(2,1,1)
    plot(x,circshift(U(:,1,j),-js),'b')
    hold on 
    plot(x,circshift(U(:,2,j),-js),'r')
    hold off
    axis([x(1) x(end) -2*a 2*a])
    grid
    set(gca,'fontsize',14)
    xlabel('x')
    title('u(x) (blue), v(x) (red)')
    subplot(2,1,2)
    uk = fft(U(:,1,j))/nx;  
    vk = fft(U(:,2,j))/nx;  
    semilogy(0:kmax,abs(uk(1:kmax+1)),'b')
    hold on
    semilogy(0:kmax,abs(vk(1:kmax+1)),'r')
    hold off
    grid
    axis([0 kmax 1e-10 1e-1])
    set(gca,'fontsize',14)
    xlabel('x')
    title('|u(k)|^2 (blue), |v(k)|^2 (red)')
    M(j) = getframe(mfig);
end
mplay(M)
gifmovie(M,'pwave_a02_k6',.05)


%uk = fft(squeeze(U(:,1,:)))/nx;
%u2 = uk(1,:);  % Eulerian mean

%usf = a^2*wp1/k1*(sin(k1*Xp-wp1*t)).^2;   % Stokes velocity  
us = a^2*wp1/(2*k1);
ueth = a^2*( wp1/(2*k1)*(cos(f*t) - 1));

figure
plot(t,ue)
hold
plot(t,ueth,'r')

% Note linear trend, implying small acceleration.  Either numerical
% issue with k=0 equations, or O(a^3) trend creeping in. Without
% that trend, it's just the IO term

dt = t(end)-t(end-1);
Xs  = us*t;  % Stokes drift
Xe  = cumsum(ue)*dt;  % Eulerian drift
Xeth = cumsum(ueth)*dt;
Xpa = mean(Xp,1);
p = polyfit(t,Xpa,1);
Xl = p(1)*t;  % estimate of actual Lagrangian drift

figure
plot(t,Xs)
hold
plot(t,Xe,'r')
plot(t,Xe+Xs,'k')
%plot(t,Xpa,'y')
plot(t,Xl,'m')
%plot(t,Xeth)
%plot(t,Xs+Xeth)
grid
set(gca,'fontsize',14)
title('Particle drift with f=1, C=1, k_0=6, a=.001, nx=256')
xlabel('t')
legend('X_s','X_e','X_s+X_e','X_l (lin trend of mean(Xp))')
%legend('X_s','X_e','X_s+X_e','X_l (lin trend of mean(Xp))','a^2\omega/(2k_0f) (sin(ft)-t)','a^2\omega/(2k_0f) (sin(ft)-t)+X_s')

%  Maybe it's not that ue is wrong .. rather, my calculation of
%  stokes is naive.  I know that higher wave modes are excited, and
%  these will add to stokes drift


% Full ue = u_2(x,t)
sig = sqrt(f^2+C^2*(2*k1)^2);
alph = wp1^3/(3*f^2*k1)+2*C^2*wp1*k1/(3*f^2)+wp1/(6*k1);
u2th = a^2*( wp1/(2*k1)*(cos(f*t) - 1) ...
     + alph*(cos(2*k1*Xp).*( cos(2*wp1*t) - cos(sig*t) )...
           + sin(2*k1*Xp).*( sin(2*wp1*t) - (2*wp1)/(sig) *sin(sig*t))));



 
% PV q = (v_x + f)/(1+h)
kf = [0:kmax -kmax-1:-1]';
q = (real(ifft(i*repmat(kf,[1 size(U,3)]).*fft(squeeze(U(:,2,:)))))+f)./(1+squeeze(U(:,3,:)));