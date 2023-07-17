% plane standing wave test

nx = 256; % must be a power of 2
%nx = 128;
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx)' - L/2;
kmax = nx/2 - 1;  % 2^n - 1
nkx = 2*kmax+1;   % -kmax:kmax
dx = L/nx;        % L_domain = 2*pi
[kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
[x_,y_] = ndgrid(x,x);

f = 1;     % = k_def if Cg = 1
Cg = 1;
aw = 0.1;  % set = 1 if v.grad v -> eps v.grad v
kw = 4;
lw = 0;
Kw2 = kw^2+lw^2;
wp = sqrt(f^2+Cg^2*Kw2);
etaw = aw*cos(kw*x_).*cos(lw*y_);
uw = aw*(kw*wp*sin(kw*x_).*cos(lw*y_)-lw*f*cos(kw*x_).*sin(lw*y_))/Kw2;
vw = aw*(lw*wp*cos(kw*x_).*sin(lw*y_)+kw*f*sin(kw*x_).*cos(lw*y_))/Kw2;
Uw(:,:,1) = uw;
Uw(:,:,2) = vw;
Uw(:,:,3) = etaw;

[U,t,ke,pe] = swk(Uw,f,Cg,10000,100);

mfig=figure('position',[700 600 700 400]);
clear M
for j=1:size(U,4)
%    xs = time(j)*c1;
%    js = floor(xs/dx);
    subplot(2,1,1)
%    plot(x,circshift(U(:,1,1,j),-js),'b')
%    hold on 
%    plot(x,circshift(U(:,1,2,j),-js),'r')
%    hold off
    plot(x,U(:,1,3,j))
    axis([x(1) x(end) -2*aw 2*aw])
    grid
    title(strcat('h(x) at t =',num2str(t(j)))) 
    set(gca,'fontsize',14)
    xlabel('x')
%    title('u(x) (blue), v(x) (red)')
    subplot(2,1,2)
%    uk = fft(U(:,1,1,j))/nx;  
%    vk = fft(U(:,1,2,j))/nx;  
    hk = fft(U(:,1,3,j))/nx;  
    semilogy(0:kmax,abs(hk(1:kmax+1)))
%    hold on
%    semilogy(0:kmax,abs(vk(1:kmax+1)),'r')
%    hold off
    grid
    axis([0 kmax 1e-20 1e-1])
    set(gca,'fontsize',14)
    xlabel('k')
%    title('|u(k)|^2 (blue), |v(k)|^2 (red)')
    title('|h(k)|^2')
    M(j) = getframe(mfig);
end
mplay(M,4)
gifmovie(M,'pwave',.05)
