% Linear Shallow water solution with stationary large-scale
% geostrophic flow, superimposed with a set of linear wave modes of
% smaller scale, with Doppler-shifting and refraction by vorticity
% of the mean flow 


% Strengths of geostrophic and wave parts
ag = 0.2;
aw = 0.1;

load ../ray_trace_sw/wavevort_231058_restart_frame100
% Contains state variable S(:,:,:) with [u,v,eta]=S(:,:,1:3)
% Use method from wavevortdecomp.m to extract geostrophic mode
% velocities ug, vg, and compute gradients spectrally
addpath ~/GoogleDrive/Matlab/rsw/
kmax = nx/2-1;
[kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
K2_   = kx_.^2 + ky_.^2;

% In .mat file, variable Cg = sqrt(g*H0).  Here we'll 
% define gH0 = Cg^2, C0 = Cg
gH0 = Cg^2;
C0 = Cg;

sig2_ = f^2 + gH0*K2_;
uk   = g2k(S(:,:,1));
vk   = g2k(S(:,:,2));
etak = g2k(S(:,:,3));
zetak  = 1i*(kx_.*vk - ky_.*uk);

etagk = (f*etak-zetak).*f./sig2_;
etag = k2g(etagk);

etag = ag*etag/max(abs(etag(:)));
etagk = g2k(etag);

zetagk = -gH0/f*etagk.*K2_;

% Extract geostrophic velocities
ugk = -1i*ky_.*(gH0/f*etagk);
vgk = 1i*kx_.*(gH0/f*etagk); 

% Get fields in x-space
etag = k2g(etagk);
ug = k2g(ugk);
vg = k2g(vgk);
vortg = k2g(zetagk);

%nx = 512; % must be a power of 2
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx)' - L/2;
[x_,y_] = ndgrid(x,x);

% Now initialize separately a spectrum of waves with random phases
k = 3:10;
l = 5:10;

% Save phases and signs of frequencies since these are constant for
% each wavemode for all times
phi = 2*pi*rand(length(k),length(l));
sign = (-1).^round(rand(length(k),length(l)));

% Time is just a parameter, but would like movie to look smooth, so
% make sure fastest wave is resolved on grid
Kmax = sqrt((max(k))^2+(max(l))^2);
Cmax = sqrt(f^2+C0^2*Kmax^2)/Kmax;
Umax = max(max(abs(ug)))+Cmax;
% choose dt = some fraction of L/Umax;

T = 50/f;
dt = .1;
t = 0:dt:T;

% For each time, create total wave field, add to geostrophic field,
% and make a movie (not saving full time dependent fields)

h = figure;

for n = 1:length(t)
    uw = 0; vw = 0; etaw = 0;
    for i = 1:length(k)
        for j = 1:length(l)
            K2 = k(i)^2+l(j)^2;
            omega = sign(i,j)*sqrt(f*(f+vortg)+C0^2*K2); 
            theta = k(i)*x_ + l(j)*y_ + phi(i,j) - (omega + k(i)*ug + l(j)*vg)*t(n);
            [u1,v1,eta1] = onewave(1,f,C0,k(i),l(j),omega,theta);
            
            uw = u1 + uw;
            vw = v1 + vw;
            etaw = eta1 + etaw;
        end
    end
    etamax = max(max(abs(etaw)));
    uw = aw*uw/etamax;
    vw = aw*vw/etamax;
    etaw = aw*etaw/etamax;

    u = ug + uw;
    v = vg + vw;
    eta = etag + etaw;
    
    % Compute vorticity
    %vort = spec2grid(-K_.^2.*grid2spec(eta));  << wrong
    vort = k2g(1i*kx_.*g2k(v)-1i*ky_.*g2k(u));
    
    clf
    %pcolor(x,x,eta'), shading interp, caxis([-aw-ag,aw+ag]), colorbar
    pcolor(x,x,vort'), shading interp, colorbar
    %set(gcf,'NextPlot','add'); set(gca,'NextPlot','add');
    %[~,hc] = contour(x,x,etag','k');
    %set(hc,'LineColor',.5*[1 1 1])
    axis image
    %title('\eta (colors) and \eta_g  (contours)')
    set(gca,'fontsize',14)
    M(n) = getframe(h);
    
end

%movie(M)


return

% Compare to output from nonlinear SW model
% Create same initial conditoin by running the above at t=0, and then:
Uin(:,:,1) = u;
Uin(:,:,2) = v;
Uin(:,:,3) = eta;

[U,t,ke,pe] = swk(Uin,f,C0,5000,100);

etabar = mean(squeeze(U(:,:,3,:)),3);
figure
contour(x,x,etabar','k'), axis image

nt = size(U,4);
h = figure;

for j=1:nt
    pcolor(x,x,U(:,:,3,j)), shading interp, axis image, caxis([-aw,aw]), colorbar;
    Ms(j)=getframe(h);
    pause(.1);
end
