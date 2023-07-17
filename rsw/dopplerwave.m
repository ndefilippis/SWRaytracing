% Linear Shallow water solution with stationary large-scale
% geostrophic flow, superimposed with a set of linear wave modes of
% smaller scale, with Doppler-shifting by the mean flow (neglecting
% refraction terms)

nx = 512; % must be a power of 2
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx)' - L/2;
[x_,y_] = ndgrid(x,x);

% kmax = nx/2 - 1;  % 2^n - 1

f = 3;   % = the deformation wavenumber k_d with C0 = 1
C0 = 1;

% Strengths of geostrophic and wave parts
ag = 0.2;
aw = 0.1;

% Geostrophic mean flow -- Childress-Soward. Parameter a controls
% how "jet-like" flow is, with a=0 being closed cells.  
a = .25;
km = 1;

etag = ag*( sin(km*x_).*sin(km*y_)+a*cos(km*x_).*cos(km*y_) );
ug   = -ag*km*C0^2/f*(sin(km*x_).*cos(km*y_)-a*cos(km*x_).*sin(km*y_));
vg   =  ag*km*C0^2/f*(cos(km*x_).*sin(km*y_)-a*sin(km*x_).*cos(km*y_));

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

T = 10/f;
dt = .1;
t = 0:dt:T;

% For each time, create total wave field, add to geostrophic field,
% and make a movie (not saving full time dependent fields)

h = figure;

ew = zeros(length(t),1);
for n = 1:length(t)
    uw = 0; vw = 0; etaw = 0;
    for i = 1:length(k)
        for j = 1:length(l)
            K2 = k(i)^2+l(j)^2;
            omega = sign(i,j)*sqrt(f^2+C0^2*K2); 
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
    
    % wave energy
    Ew = uw.^2 + vw.^2 + C0^2*etaw.^2;
    ew(n) = sum(sum(Ew));
    
    u = ug + uw;
    v = vg + vw;
    eta = etag + etaw;
    
    clf
    %pcolor(x,x,eta'), shading interp, caxis([-aw-ag,aw+ag]), colorbar
    pcolor(x,x,Ew'), shading interp, colorbar
    %set(gcf,'NextPlot','add'); set(gca,'NextPlot','add');
    %[~,hc] = contour(x,x,etag','k');
    %set(hc,'LineColor',.5*[1 1 1])
    axis image
    %title('\eta (colors) and \eta_g  (contours)')
    set(gca,'fontsize',14)
    pause(0.01)
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
