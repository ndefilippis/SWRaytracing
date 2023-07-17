function [zeta,q,B] = get_swvort(u,v,h,hb,dx,dy,Cgsq,Roi,Be)

%  [zeta,q,B] = get_swvort(u,v,h,hb,dx,dy,Cgsq,Roi,Be)
%  Get SW vorticity zeta, PV q and Bernouli B from u, v and h. hb
%  is topography field (set to 0 if none).  Cgsq = g*H_0/U^2, Roi =
%  f_0*L/U and Be = beta_0*L^2/U.

ny = size(v,2);

H = h-hb;

% Coriolis parameter
if (Be==0)
  f = Roi;
else
  y = linspace(-1/2,1/2,ny);
  fvec = Roi + Be*y;
  f =  meshgrid(fvec,1:nx-1);           % f on v grid
end

% Get vorticity (i,j) -- I think below is free-slip - no vorticity
% at boundaries
zeta = zeros(size(u,1),size(v,2));
zeta(2:end-1,2:end-1) = diff(v(:,2:end-1),1,1)/dx - diff(u(2:end-1,:),1,2)/dy;

% Get Bernouli on h points
B = Cgsq*h + ((avg(u,1)).^2 + (avg(v,2)).^2)/2;

% Get PV on zeta points
q = zeros(size(zeta));
q = (f(:,2:end-1) + zeta(2:end-1,2:end-1))./H;
