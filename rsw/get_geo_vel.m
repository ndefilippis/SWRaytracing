function [u,v] = get_geo_vel(h,dx,dy,f);

% [u,v] = get_geo_vel(h,dx,dy,f)
% Get SW geostrophic velocities from height field h.  Returns u and
% v on a C-grid relative to h, hence size(u) = [nx-1,ny] and
% size(v) = [nx,ny-1].  f is coriolis f on vector congruent with
% second dimension of h (y dim)..

nx = size(h,1)+1;
ny = size(h,2)+1;

fh = (f(2:ny)+f(1:ny-1))/2;  % f on staggered grid, for
                               % calculation of v

f_u  = meshgrid(avg(f,2),1:nx);      % f on u grid
f_v =  meshgrid(f,1:nx-1);           % f on v grid

dhx = zeros(size(h)); dhy = zeros(size(h));

dhy(:,2:ny-1) = (h(:,3:ny)-h(:,1:ny-2))/(2*dy);  % Get dh/dy on h points (i,j)
dhy(:,1) = (h(:,2)-h(:,1))/dy;
dhy(:,ny) = (h(:,ny)-h(:,ny-1))/dy;

u = -(dhy(2:nx,:)+dhy(1:nx-1,:))./(2*f_u(1:nx-1,:)); % Avg dh/dy onto u pts
                                                


dhx(2:nx-1,:) = (h(3:nx,:)-h(1:nx-2,:))/(2*dx);  % Get dh/dx on h points 
dhx(1,:) = (h(2,:)-h(1,:))/dx;
dhx(nx,:) = (h(nx,:)-h(nx-1,:))/dx;

v = (dhx(:,2:ny)+dhx(:,1:ny-1))./(2*f_v);       % Avg dh/dx onto v
                                                 % points (i,j+1/2)
