% Notes:  a run of length 1/(f*Fr^2) is very short compared to the
% length of the nonlinear SW run on the HPC.  With f=3, end time of
% run-2371058 is t_end*f = 2328 (and dt=7.7614e-04).  To do:
% * Check nondim scaling of swk.m
% * Figure out how quickly wavesignal evolves, using shorter time
% frequency/wavenumber spectrum.  Possibly need to run model till
% reasonable geostrophic field, kill waves, then restart. I may
% have done this already (check).

% Use frame from wavevort simulations on HPC to advect packets
load wavevort_231058_restart_frame100
% Contains state variable S(:,:,:) with [u,v,eta]=S(:,:,1:3)
% Use method from wavevortdecomp.m to extract geostrophic mode
% velocities ug, vg, and compute gradients spectrally
addpath ../rsw/
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
zetagk = -gH0/f*etagk.*K2_;

% Extract geostrophic velocities
ugk = -1i*ky_.*(gH0/f*etagk);
vgk = 1i*kx_.*(gH0/f*etagk); 

% Get gradients
ugxk = 1i*kx_.*ugk;
ugyk = 1i*ky_.*ugk;
vgxk = 1i*kx_.*vgk;
vgyk = 1i*ky_.*vgk;

% Get fields in x-space
etag = k2g(etagk);
H = 1+etag;
U.u = k2g(ugk);
U.v = k2g(vgk);

GradU.u_x = k2g(ugxk);
GradU.u_y = k2g(ugyk);
GradU.v_x = k2g(vgxk);
GradU.v_y = k2g(vgyk);
clear kx_ ky_ K2_ sig2_ uk vk etak zetak ugk vgk ugxk ugyk vgxk vgyk

figure
pcolor(GradU.u_y), shading interp, axis image
% Unfortunately gradients have a lot of finescale striations
% ... perhaps smooth them?  First we'll ignore them and see what happens


% Added Aug 5, 21:  
% For RSW with mean H(x,y) in geostrophic balance with [U,V](x,y), so
% that div(U) = 0, have 
% omega = sqrt(f^2 + g*H(x,y)*K^2)
% and
% C = [c_x,c_y] = g*H(x,y)/omega(k,l,x,y) * [k,l]
% Include equation for wave action:
% a = |A|^2 omega as
% da/dt = -a*div_x(C)
% Compute div_x(C) at start, since it's fixed
% div_x(C) = omega^(-1) * ( k*g*H_x + l*g*H_y - c_x^2 - c_y^2)
% From the SW output above:  H = etag, [U,V] = [U.u,U.v]
% and from this analysis, note that everywhere I've written g it
% should really be g*H0 = Cg^2 since H = etag is nondimensionalized
% by H0.  But continuing with notation above, g*[H_x,H_y] = f*[V,-U]   
% so
% div_x(C) = omega^(-1) * ( k*f*V - l*f*U - c_x^2 - c_y^2)
% thus no need to compute any more derivatives


% Parameters
%f = 3;          set by load command
%nx = 512;       

L = 2*pi;       % Domain scale (so wavenumbers are integers)
kd = f/C0;    % Deformation wavenumber (= f with C0 = 1) 
ki = kd*10;     % Initial packet wavenumber magnitude
omega0 = sqrt( f^2 + gH0*ki^2 );% = 10*f with f=3, gH0=1, ki=kd*10 

% Estimate Fr
speed = sqrt(U.u.^2+U.v.^2);
figure
pcolor(speed/C0), shading interp, axis image, colorbar
% max = .25, cool
% Set U0 as maximum speed for CFL
U0 = max(speed(:));
Ro = U0*kd/f   % Rossby # (= Fr with L=1/kd)
Fr = U0/C0
np = 10;       % Number of packets, changed from 100 to 5, 2023/04/10

dx = L/nx;      % Grid resolution     
dt = .3*dx/max(C0,U0);  % Time resolution

% 2023/04/10: Changed to 0.1 
Tend = 20/(f*Fr^2); 
nsteps = round(Tend/dt)  % Number of timesteps for simulation

% Time and space grids
t = 0:dt:dt*(nsteps-1);
x = 0:dx:dx*(nx-1);
[x_,y_] = ndgrid(x,x);

% Initial packet positions and wavenumbers
clear P
P(np,nsteps) = struct();
for i = 1:np
    P(i,1).x = rand*L;
    P(i,1).y = rand*L;
    P(i,1).k = ki*cos(2*pi*i/np);
    P(i,1).l = ki*sin(2*pi*i/np);
    P(i,1).a = 1;
end    

% Advect packets
for i = 1:np
    for j = 2:nsteps
        %P(i,j) = step_packet(P(i,j-1),U,GradU,sqrt(gH0),f,dx,dx,dt);
        P(i,j) = step_packet_xka(P(i,j-1),U,GradU,H,C0,f,dx,dx,dt);
    end
end

% Save group velocities, frequency, wavenumber for each particle
clear omega C Cmag K
for i=1:np
    for j=1:nsteps
        [C(i,j),omega(i,j)] = cg_sw(P(i,j).k,P(i,j).l,C0,f,H);
        Cmag(i,j) = sqrt(C(i,j).x^2+C(i,j).y^2);
        K(i,j) = sqrt(P(i,j).k^2+P(i,j).l^2);
    end
end

%[C,omega,divC,gradomega] = cg_sw(P.k,P.l,C0,f,H,U);

%% Plots 

figure
semilogy(t/f,mean(omega,1)/f);
grid

% Check dispersion relation
%figure
%plot(sqrt(mean(K,1).^2/kd^2),sqrt(mean(omega,1).^2/f^2-1))

% Plot evolution of rays in x,y space, overlayed on streamfunction
figure
% Below is a hack to make a controuf plot partially transparent
[~, hContour] = contourf(x,x,etag'); 
axis image
drawnow;  
hFills = hContour.FacePrims;  % array of TriangleStrip objects
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
for idx = 1 : numel(hFills)
   hFills(idx).ColorData(4) = 50;   % default=255
end
% Now add trajectories on top of streamfunction
hold on
plot([P(:,1).x],[P(:,1).y],'g.','Markersize',30)
for j=2:5:nsteps-1
    plot(mod([P(:,j).x],L),mod([P(:,j).y],L),'k.','Markersize',.5)
end
plot(mod([P(:,end).x],L),mod([P(:,end).y],L),'r.','Markersize',30)

% Plot evolution of rays in k,l space
kvec = 1:max(K(:))/kd;
figure
plot([P(:,1).k]/kd,[P(:,1).l]/kd,'g.','Markersize',30)
axis image
hold on
for j=2:2:nsteps-1
    plot([P(:,j).k]/kd,[P(:,j).l]/kd,'k.','Markersize',.5)
end
plot([P(:,end).k]/kd,[P(:,end).l]/kd,'r.','Markersize',30)

% Plot trajectory, with velocity vectors U+C at each point
figure
for ip = 1:10:np
   plot(mod([P(ip,:).x],L),mod([P(ip,:).y],L),'k.','Markersize',.5)
   hold on
   for j=1:5:nsteps
       u = interpolate(P(ip,j).x,P(ip,j).y,U.u,dx,dx);
       v = interpolate(P(ip,j).x,P(ip,j).y,U.v,dx,dx);
       quiver([P(ip,j).x],[P(ip,j).y],[C(ip,j).x+u],[C(ip,j).y+v],.2,'b')
        quiver([P(ip,j).x],[P(ip,j).y],u,v,.5,'b')   
   end
end

% Observations
%   * anisotropy (via parameter a) greatly increases frequency
%     diffusion, even at low Fr
%   * 

% Plot simulated wavefield using fiduciary form for A(x)

figure

% for j=1:nsteps
%     etaw = zeros(size(H));
%     for i=1:np
%         Hl = H(ceil(mod(P(i,j).x/dx, nx)),ceil(mod(P(i,j).y/dx, nx)));
%         [~,w] = cg_sw(P(i,j).k,P(i,j).l,C0,f,Hl);        
%         Amax = sqrt(P(i,j).a/w);
%         A = ampfunc(x_,y_,P(i,j).x,P(i,j).y,Amax,2*pi/50);
%         etat = A.*cos(P(i,j).k*x_ + P(i,j).l*y_ - w*t(j));
%         etaw = etaw + etat;
%     end
%     pcolor((etag+etaw)'), shading interp, axis image, colorbar, caxis([-.5,.5])
%     pause(.01)
% end
