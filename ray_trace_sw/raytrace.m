
% Parameters
C0 = 1;         % GW wave speed (sqrt(gH))
Fr = .1;         % Froude # = U0/C0
U0 = Fr;        % Mean flow speed
f = 4;          % Nondim Coriolis
L = 2*pi;       % Domain scale (so wavenumbers are integers)
kd = f/C0^2;    % Deformation wavenumber (= f with C0 = 1) 
km = kd;        % Mean flow wavenumber
Ro = U0*kd/f;   % Rossby # (= Fr with L=1/kd)
a = .25;        % Shear parameter for Childress-Soward flow
%a = 0;
ki = km*10;     % Initial packet wavenumber magnitude
omega0 = sqrt( f^2 + C0^2*ki^2 );% = 10*f with f=3, C0=1, ki=kd*10 

np = 5;       % Number of packets
nx = 256;       % Spatial grid size

dx = L/nx;      % Grid resolution     
dt = .3*dx/max(C0,U0);  % Time resolution

Tend = 1/(f*Fr^2);
nsteps = round(Tend/dt);  % Number of timesteps for simulation

% Time and space grids
t = 0:dt:dt*(nsteps-1);
x = 0:dx:dx*(nx-1);
[x_,y_] = ndgrid(x,x);

% Childress-Soward flow
psi = U0*km^(-1)*( sin(km*x_).*sin(km*y_)+a*cos(km*x_).*cos(km*y_) );
U.u = -U0*(sin(km*x_).*cos(km*y_)-a*cos(km*x_).*sin(km*y_));
U.v =  U0*(cos(km*x_).*sin(km*y_)-a*sin(km*x_).*cos(km*y_));
GradU.u_x = -km*U0*(cos(km*x_).*cos(km*y_)+a*sin(km*x_).*sin(km*y_));
GradU.u_y = km*U0*(sin(km*x_).*sin(km*y_)+a*cos(km*x_).*cos(km*y_));
GradU.v_x = -km*U0*(sin(km*x_).*sin(km*y_)+a*cos(km*x_)*cos(km*y_));
GradU.v_y = km*U0*(cos(km*x_).*cos(km*y_)+a*sin(km*x_).*sin(km*y_));
vort = GradU.v_x - GradU.u_y;

% Initial packet positions and wavenumbers
clear P
P(np,nsteps) = struct();
for i = 1:np
    P(i,1).x = rand*L;
    P(i,1).y = rand*L;
    P(i,1).k = ki*cos(2*pi*i/np);
    P(i,1).l = ki*sin(2*pi*i/np);
end    

% Advect packets
for i = 1:np
    for j = 2:nsteps
        P(i,j) = step_packet(P(i,j-1),U,GradU,C0,f,dx,dx,dt);
    end
end

% Save group velocities, frequency, wavenumber for each particle 
for i=1:np
    for j=1:nsteps
        [C(i,j),omega(i,j),omega_abs(i,j)] = cg_sw(P(i,j).k,P(i,j).l,C0,f,U);
        Cmag(i,j) = sqrt(C(i,j).x^2+C(i,j).y^2);
        K(i,j) = sqrt(P(i,j).k^2+P(i,j).l^2);
    end
end

%% Plots 

figure
semilogy(t/f,mean(omega,1)/f);
grid

% Check dispersion relation
figure
plot(sqrt(mean(K,1).^2/kd^2),sqrt(mean(omega,1).^2/f^2-1))

% Plot evolution of rays in x,y space, overlayed on streamfunction
figure
% Below is a hack to make a controuf plot partially transparent
[~, hContour] = contourf(x,x,psi'); 
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
%figure
%for ip = 1:10:np
%    plot(mod([P(ip,:).x],L),mod([P(ip,:).y],L),'k.','Markersize',.5)
%    hold on
%    for j=1:5:nsteps
%        u = interpolate(P(ip,j).x,P(ip,j).y,U.u,dx,dx);
%        v = interpolate(P(ip,j).x,P(ip,j).y,U.v,dx,dx);
%        quiver([P(ip,j).x],[P(ip,j).y],[C(ip,j).x+u],[C(ip,j).y+v],.2,'b')
        %quiver([P(ip,j).x],[P(ip,j).y],u,v,.5,'b')   
%    end
%end

% Observations
%   * anisotropy (via parameter a) greatly increases frequency
%     diffusion, even at low Fr
%   * 