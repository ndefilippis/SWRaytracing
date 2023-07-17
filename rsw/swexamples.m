
nx = 256; % must be a power of 2
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx) - L/2;
[x_,y_] = ndgrid(x,x);
kmax = nx/2 - 1;  % 2^n - 1
nkx = 2*kmax+1;   % -kmax:kmax
nky = kmax+1;     % 0:kmax
dx = 2*pi/nx;     % L_domain = 2*pi
[ikx_,iky_] = ndgrid(-kmax:kmax,0:kmax);
K_   = sqrt(ikx_.^2 + iky_.^2);
ikx_ = sqrt(-1)*ikx_; 
iky_ = sqrt(-1)*iky_;

% Make a wave packet, oriented by theta from x:

f = 1;
Cg = 1;
h0 = .1;
k0 = 5;
theta = 0;

w = -sqrt(f^2+Cg^2*k0^2);

xp_ = x_*cos(theta)-y_*sin(theta);
yp_ = x_*sin(theta)+y_*cos(theta);
xp0 = pi/6;
yp0 = -pi/4;
henv = exp(-(10*(xp_-xp0)/L).^2-(10*(yp_-yp0)/L).^2);
%henv=1;
hin = h0*henv.*cos(k0*xp_);

Sin(:,:,1)=h0*w/k0*henv.*cos(k0*x_);
Sin(:,:,2)=h0*f/k0*henv.*sin(k0*x_);
Sin(:,:,3)=hin;

[Sout,t,ke,pe,Xp,hmov] = swk(Sin,f,Cg,5000,10);

%movie(hmov,1,3)

%-------------------------------------------------------------
% 0 PV geostrophic adjustment

f = 10;
Cg = 10;
h0 = .01;
b = 10;

% Make a localized jump in surface, and set vorticity = f*h, so
% that PV = (f+vort)/(1+h) = f

henv = exp(-(b*(y_)/L).^2);
hin = h0*henv.*x_./(x_.^4+.01);
psik = -g2k(hin)./K_.^2; 
psik(kmax+1,1)=0;
u = f*k2g(-iky_.*psik);
v = f*k2g(ikx_.*psik); 

Sin(:,:,1)=u;
Sin(:,:,2)=v;
Sin(:,:,3)=hin;

np = 20;
Xpin.x = .2*randn(np,1);
Xpin.y = .2*randn(np,1);

[Sout,t,ke,pe,Xp,hmov] = swk(Sin,f,Cg,1000,25,Xpin);

% blows up at frame 40, but PV at 20 is still constant:

zeta = k2g(ikx_.*g2k(Sout(:,:,2,20))-iky_.*g2k(Sout(:,:,1,20)));
q=(f+zeta)./(1+Sout(:,:,3,20));
figure
plot(q)  % plot all rows of q 
axis([0 140 9 11])




%-------------------------------------------------------------
% For geostrophic IC: 

Cg=10;
f=50;
r=Cg^2/f;

hin = exp(-(10*x_/L).^2-(20*y_/L).^2); %.*cos(40*x_);
hin = hin-mean(hin(:));
hk = g2k(hin);
ug = -r*k2g(iky_.*hk);
vg =  r*k2g(ikx_.*hk);
   
% Normalize to make max(abs(velocity))=1
U0 = max([max(abs(ug(:))) max(abs(vg(:)))]);
ug=ug/U0; vg=vg/U0; hin=hin/U0;

Sin(:,:,1)=ug;
Sin(:,:,2)=vg;
Sin(:,:,3)=hin;

np = 20;
Xpin.x = .2*randn(np,1);
Xpin.y = .2*randn(np,1);

[Sout,t,ke,pe,Xp,hmov] = swk(Sin,f,Cg,5000,25,Xpin);

xp=[Xp.x]; yp=[Xp.y];
j=100;
divu = k2g(ikx_.*g2k(Sout(:,:,1,j))+iky_.*g2k(Sout(:,:,2,j)));
plotfield(divu,1,'pc')
uadv = ug.*k2g(ikx_.*g2k(ug))+vg.*k2g(iky_.*g2k(ug));
plotfield(uadv,1,'pc')
plotfield(f*ug,1,'pc')
% nonlinear terms are largish here

%-------------------------------------------------------------
% Vortex in rigid lid limit
Cg=0;
f=0;
hin = exp(-(6*x_/L).^2-(24*y_/L).^2); % vorticity
hk = -(1./K_).^2.*g2k(hin);           % psi
hk(kmax+1,1) = 0;  % remove k,l = 0,0 compt
ug = -k2g(iky_.*hk);
vg =  k2g(ikx_.*hk);

Sin(:,:,1)=ug;
Sin(:,:,2)=vg;
Sin(:,:,3)=0;

np = 20;
Xpin.x = .2*randn(np,1);
Xpin.y = .2*randn(np,1);
% manually alter plotstuff in code to plot vorticity instead of h

[Sout,t,ke,pe,Xp,hmov] = swk(Sin,f,Cg,1000,25,Xpin);

%-------------------------------------------------------------
% I.O.

% v = cos(ft) => u = sin(ft);
f=1;

Sin=zeros(nx,nx,3);
Sin(:,:,2)=1;
[Sout,t,ke,pe,Xp] = swk(Sin,f,Cg,1000,25,Xpin);

figure
plot(t*f/(2*pi),squeeze(Sout(1,1,1,:)))
hold
plot(t*f/(2*pi),squeeze(Sout(1,1,2,:)))
plot(t*f/(2*pi),squeeze(Sout(1,1,3,:)))
grid
legend('u','v','h')

%-------------------------------------------------------------
% Plane gravity wave

f=1;
Cg=1;
h0=.1;
k0=4;
wp = sqrt(f^2+Cg^2*k0^2);

Sin(:,:,1)=h0*wp/k0*cos(k0*x_);
Sin(:,:,2)=h0*f/k0*sin(k0*x_);
Sin(:,:,3)=h0*cos(k0*x_);

% 7/21/17:  removed np and nutune from arg list;  added ability to
% read in list of initial particle positions

np = 20;
Xpin.x = .2*randn(np,1);
Xpin.y = .2*randn(np,1);

[Sout,t,ke,pe,Xp,hmov] = swk(Sin,f,Cg,2000,25,Xpin);

plot(Xp.x(1,:),Xp.y(1,:))
clear M
for j=1:size(Sout,4)
    plot(x,squeeze(Sout(:,1,3,j)))
    grid
    %axis([-pi pi -.1 .1])
    axis tight
    drawnow 
    M(j)=getframe;
end

% Stokes drift:  
us = h0^2*wp/(2*k0);

%-------------------------------------------------------------
% Plane gravity wave

f=0;
Cg=1;
h0=.1;
k0=2;
wp = Cg*k0;

Sin(:,:,1)=h0*wp/k0*cos(k0*x_);
Sin(:,:,2)=h0*f/k0*sin(k0*x_);
Sin(:,:,3)=h0*cos(k0*x_);

[Sout,t,ke,pe,hmov] = swk(Sin,f,Cg,5000,25);

for j=1:size(Sout,4)
    plot(x,squeeze(Sout(:,1,3,j)))
    grid
    %axis([-pi pi -.1 .1])
    axis tight
    drawnow 
    M(j)=getframe;
end

movie(M,1,3)

%-------------------------------------------------------------
% For tsb

Ro = .0125;
Bu = 10*Ro;
h0=.1;
k0=2;
f = 1;
wp = sqrt(f^2+Cg^2*k0^2);

Sin = zeros(nx,nx,3);
Sin(:,:,1)=h0*wp/k0*cos(k0*x_);
Sin(:,:,2)=h0*f/k0*sin(k0*x_);
Sin(:,:,3)=0;

Psi = sin(x_).*sin(y_);
[Sout,time,hmov] = tsbn(Sin,Bu,Ro,Psi,500,25);

for j=1:size(Sout,4)
    plot(x,squeeze(Sout(:,1,3,j)))
    grid
    %axis([-pi pi -.1 .1])
    axis tight
    drawnow 
    M(j)=getframe;
end

%----------------------------------------------------------------
% 2 Plane gravity waves with opposite phase velocity

f=1;
Cg=1;
% wp = sqrt(f^2+Cg^2*k0^2);
nx = 256; % must be a power of 2
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx);
[x_,y_] = ndgrid(x,x);

h1=.05; h2=.1;
% k1X=5; k1Y=7;  k2X=7; k2Y=5; 
k1X=1; k1Y=0;  k2X=1; k2Y=0; 

K1=sqrt(k1X^2+k1Y^2); K2=sqrt(k2X^2+k2Y^2); 

omega1=sqrt(1+K1^2); omega2=sqrt(1+K2^2); 

h0= h1*exp(1i.*(k1X.*x_+k1Y.*y_)) + h2*exp(1i.*(k2X.*x_+k2Y.*y_)) ;
u0= (omega1*k1X + 1i*k1Y)*(h1/K1^2).*exp(1i.*(k1X.*x_+k1Y.*y_)) + (-omega2*k2X + 1i*k2Y)*(h2/K2^2).*exp(1i.*(k2X.*x_+k2Y.*y_)) ;
v0= (omega1*k1Y - 1i*k1X)*(h1/K1^2).*exp(1i.*(k1X.*x_+k1Y.*y_)) + (-omega2*k2Y - 1i*k2X)*(h2/K2^2).*exp(1i.*(k2X.*x_+k2Y.*y_)) ;

h0=h0 + conj(h0); u0=u0 + conj(u0); v0=v0 + conj(v0);
 
Sin(:,:,1)=u0;     
Sin(:,:,2)=v0;      
Sin(:,:,3)=h0;     

[Sout,t,ke,pe,hmov] = swk(Sin,f,Cg,1000,25);

%----------------------------------------------------------------
% Plane gravity wave -- Stokes drift

f  = 1;
Cg = 2;
h0 = .01;
k0 = 2;  % initial wave in x direction

wp = sqrt(f^2+Cg^2*k0^2);  % 

% With l0 = 0, expect Stokes trajectory
% [Xs,Ys] = -h0/k0 * [sin(k0*x-wp*t),-f/wp * cos(k0*x-wp*t)]  
% and drift is (averaging over wave periods)
% Us = wp*h0^2/(2*k0), Vs = 0

Sin(:,:,1) = h0*wp/k0*cos(k0*x_);
Sin(:,:,2) = h0*f/k0*sin(k0*x_);
Sin(:,:,3) = h0*cos(k0*x_);

[Sout,t,ke,pe,hmov,Xp] = swk(Sin,f,Cg,3000,25,1,5);

plot(Xp.x(1,:),Xp.y(1,:))

for j=1:size(Sout,4)
    plot(x,squeeze(Sout(:,1,3,j)))
    grid
    %axis([-pi pi -.1 .1])
    axis tight
    drawnow 
    M(j)=getframe;
end

ts = runavg(t,5);
Xs = runavg(Xp.x(1,:),5);
Us = Cg*h0^2/2;

figure
plot(t,Xp.x(1,:));
hold
plot(ts,Xs,'r')
grid
plot(t,Us*t,'k')
set(gca,'fontsize',16)
xlabel('t')
title('Trajectories')
legend('X(t)','<X>(t)','X_s(t) = (0.5) C_g h0^2','location','southeast')