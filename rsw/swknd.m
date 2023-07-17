function [Sout,time,ke,pe,hmov,Xp] = swknd(Sin,ep,gam,numsteps,savestep,nutune,np)

%  [Sout,time,ke,pe,hmov,Xp] = swk(Sin,ep,gam,numsteps,savestep,nutune,np)
%
%  Solves rotating shallow water equations with H = gam*(1+ep*h)
%
%  u_t - v*(1 + ep*zeta) + B_x = nu*del^(a) u
%  v_t + u*(1 + ep*zeta) + B_y = nu*del^(a) v
%  h_t + gam*[(1+ep*h)*u]_x + gam*[(1+ep*h)*v]_y = 0 
%  
%  zeta = v_x-u_y 
%  B = gam*[ep*(u^2+v^2)/2 + h]    
%
%  Inputs
%
%  Sin:       N x N x 3 array containing initial u, v, h, respectively
%  ep:        U/f*Ld, where Ld = sqrt(g*H_0)/f
%  gam:       Ld/L = sqrt(Bu)
%  numsteps:  Total number of timesteps 
%  savestep:  Frequency, in timesteps, to save output    
%  np:        Advect np^2 particles (optional)
%
%  Outputs
%
%  Sout:      Arranged as Sin, but with 4th dimension for time
%  time:      Times at which output is saved
%  ke:        Time series of KE
%  pe:        Time series of PE
%  hmov:      Movie of h field (if hmov included in output list)
%  Xp:        Structure with coordinates Xp.x, Xp.y of particles
%    
%  Numerical details
%
%  Model is spectal, in square domain of size 2*pi x 2*pi.  Input
%  fields must have N = 2^n, where n is an integer.  Nonlinear terms
%  are done in physical space using dealiased product via Orszag
%  method.  Uses AB3 timestepping with trapezoidal hyperviscosity of
%  order a.  Timestep and hyperviscosity are set adaptively, via 
%  dt = dttune*dx/max(|u|) and nu = nutune*dx^a*zeta_rms.
%
%  Tuning factors dttune and nutune and hyperviscous order a
%  can be set by editing this file directly.

% Tuning
a = 8;         % hyperviscosity order, nu del^a u 
dttune = .5;   % Courant number
%nutune = 2;

% Advect particles?
if (nargin>6), particles=true; end

% Check for outputs requested
makemov=false;
if (nargout>4), makemov=true;  end
   
% Set global params for use in rhs functions
global nx ikx_ iky_ kmax nkx nky u v h zeta divuk 
global kgfac gkfac1 gkfac2 damask 
global xp yp t

% Maximum phase speed.  Note omega = +/-sqrt(eps*gam^2*K^2+1) 
Cmax = sqrt(ep*gam^2+1);

% Get and check dimensions
[nx,ny] = size(Sin(:,:,1));
if (nx~=ny), error('must have nx = ny'); end
if (mod(log2(nx),1)~=0), error('must have nx = 2^n, n integer'); end

% Dimensions of spectral fields -- stored on upper half-plane in kx,ky
% (lower half-plane given by conjugate symmetry).
kmax = nx/2 - 1;  % 2^n - 1
nkx = 2*kmax+1;   % -kmax:kmax
nky = kmax+1;     % 0:kmax

L = 2*pi;
dx = L/nx;

% Set up arrays of wavenumbers for computing derivatives
[ikx_,iky_] = ndgrid(-kmax:kmax,0:kmax);
K_   = sqrt(ikx_.^2 + iky_.^2);
ikx_ = sqrt(-1)*ikx_; 
iky_ = sqrt(-1)*iky_;

% Initialize fields for spectral products with Orszag dealiasing
kcut      = sqrt(8./9.)*(kmax+1);
damask    = ones(size(ikx_));
damask(K_>kcut)    = 0.;
damask(1:kmax,1)   = 0.;  % changed from 1:kmax+1 to keep 0,0 value.

eipik  = exp(pi*(ikx_+iky_)/nx);
kgfac  = 1 + sqrt(-1)*fulspec(eipik);
gkfac1 = (1 - sqrt(-1)*conj(eipik))/4;
gkfac2 = (1 + sqrt(-1)*conj(eipik))/4;

% For hyperviscous filter
Ka_ = K_.^a;  
clear K_ eipik 

% For particle advection:  number of particles = np^2
% Initialize on uniform grid
if (particles)  
    x0 = (0:np-1)/np * L + 1e-7;
    [foox,fooy] = ndgrid(x0,x0);  % positions of np^2 paticles at t=0
    xp = foox(:);  yp = fooy(:);  % use 1-d arrays
end

% Get initial spectral PV
Sk = g2k(Sin);

% Set params for AB3+trapezoidal diffusion (Durran 3.81)
a1 = 23/12;   a2 = -16/12;  a3 = 5/12;  

% Preallocate arrays to hold saved output
nframes = floor(numsteps/savestep);
Sout = zeros(nx,ny,3,nframes);  % output on grid
time = zeros(1,nframes); 
ke = zeros(1,nframes); 
pe = zeros(1,nframes); 
if (particles), Xp.x = zeros(np^2,nframes); Xp.y=Xp.x; end

% Set counters and temp fields
frame = 0;  t = 0;  n = 0; 
Rk = 0; Rkm1 = 0; Rkm2 = 0;

keepgoing = true;
while keepgoing
        
    % Save n-1 and n-2 rhs and get next one.  
    % getrhs sets u, v, h, too, but they have imaginary parts
    % holding staggered grid fields; use real() for computations etc 
    Rkm2 = Rkm1;
    Rkm1 = Rk;
    Rk = getrhs(Sk,ep,gam);  
    
    if (n==0) Rkm1 = Rk; Rkm2 = Rk; end

   %     Umax = max([max(max(abs(real(u)))) max(max(abs(real(v)))) Cmax]);
   Umax = max([max(abs(real(u(:)))) max(abs(real(v(:)))) Cmax]);

    % Exit if blowing up, and save last field.
    if (Umax>1e6)  
        disp(strcat('Blow up! Umax= ',num2str(Umax),', t= ',num2str(t)))
        Sout(:,:,:,frame+1) = k2g(Sk); 
        keepgoing = false;
    end
    
    % Adapt dt and nu
    dt = dttune*dx/Umax;   % Courant condition 
    nu = nutune*dx^a/dt;

    % Trapezoidal diffusion operators
    filterup = (1 - (dt/2)*nu*Ka_);                
    filterdn = 1./(1 + (dt/2)*nu*Ka_);
    
    % Save output at frequency savestep
    if (mod(n,savestep)==0||n==0)  
        frame = frame+1;  
        ke(frame) = .5*sum(real(1+ep*h(:)).*(real(u(:)).^2+real(v(:)).^2));
        pe(frame) = .5/ep^2*sum(real(1+ep*h(:)).^2);
        %hbar(frame) = sum(real(h(:)))/nx^2;
        Sout(:,:,1,frame) = real(u);    
        Sout(:,:,2,frame) = real(v);    
        Sout(:,:,3,frame) = real(h);  
        time(frame) = t;
        disp(strcat('Wrote frame >',num2str(frame),' out of >',num2str(nframes)))
        disp(strcat('max(|u|) = ',num2str(Umax),', dt = ',num2str(dt),', nu = ',num2str(nu)))
        if (makemov)
            hc = plotstuff(frame,particles,ep,gam);
            hmov(frame) = getframe(hc);
        end
        if (particles)
            Xp.x(:,frame)=xp;
            Xp.y(:,frame)=yp;
        end
    end
    
    % Timestep and diffuse
    Sk = repmat(filterdn,[1 1 3]).*(repmat(filterup,[1 1 3]).*Sk + dt*(a1*Rk + a2*Rkm1 + a3*Rkm2));

    if (particles) % advect particles
        [xp,yp] = advect_particles(xp,yp,real(u),real(v),dt,dx);
    end
    
    n = n+1;
    t = t+dt;  % clock
    
    if (n==numsteps), disp('End reached'), keepgoing=false; end

end

return

%-------------------------------------------------------------------
% Internal functions
%-------------------------------------------------------------------

function Rk = getrhs(Sk,ep,gam)
   
    global ikx_ iky_ u v h zeta divuk

    u = k2gp(Sk(:,:,1));
    v = k2gp(Sk(:,:,2));
    h = k2gp(Sk(:,:,3));

    zeta = k2gp(ikx_.*Sk(:,:,2) - iky_.*Sk(:,:,1));     % vorticity
    divuk = ikx_.*Sk(:,:,1) + iky_.*Sk(:,:,2);
    Bk = 0.5*gam*ep*gp2k(gprod(u,u)+gprod(v,v)) + gam*Sk(:,:,3);  % Bernouli in k-spc
    
    Rk(:,:,1) =  gam*ep*gp2k(gprod(v,zeta)) + Sk(:,:,2) - ikx_.*Bk;
    Rk(:,:,2) = -gam*ep*gp2k(gprod(u,zeta)) - Sk(:,:,1) - iky_.*Bk;
    Rk(:,:,3) = -gam*( ep*ikx_.*gp2k(gprod(u,h)) +ep*iky_.*gp2k(gprod(v,h)) +divuk);
    
    return
    
%-------------------------------------------------------------------
        
function fg = k2gp(fk)

   % Transform to grid space, packing fields shifted by dx/2 into
   % imaginary part
    
   global nx kgfac damask

   fkt = fulspec(damask.*fk).*kgfac;
   fg  = nx^2*ifft2(ifftshift(fkt));

   return
   
%-------------------------------------------------------------------

function prodg = gprod(f,g)
    
    % Product of real parts and imaginary parts, loaded
    % respectively into real and imaginary part of output
        
    prodg = real(f).*real(g) + sqrt(-1)*imag(f).*imag(g);

    return
    
%-------------------------------------------------------------------

function prodk = gp2k(prodg)
    
    % Transform grid field, holding product of fields, 
    % with imaginary part holding field
    % shifted by dx/2, to k-space, and average result
        
    global nx kmax gkfac1 gkfac2
  
    Wk  = fftshift(fft2(prodg))/nx^2;
    
    % Extract spectral products on grid and shifted grid, and average.
    
    Wk_up = Wk(2:end,kmax+2:end);
    Wk_dn = rot90(rot90((conj(Wk(2:end,2:kmax+2)))));
    %prodk = ((1 - icalphak).*Wk_up + (1 + icalphak).*Wk_dn)/4;
    prodk = gkfac1.*Wk_up + gkfac2.*Wk_dn;

    return

%-------------------------------------------------------------------
    
function fk = g2k(fg)
    
    % Just for transforming gridded input field
        
    global nx kmax
    
    for j=1:size(fg,3)
        fkt = fftshift(fft2(fg(:,:,j)))/nx^2;
        fk(:,:,j) = fkt(2:end,kmax+2:end);      
    end
    
    return
    
%-------------------------------------------------------------------
        
function fg = k2g(fk)

   % Transform to grid space - just for output.
    
   global nx 

   for j=1:size(fk,3)
       fg(:,:,j)  = nx^2*ifft2(ifftshift(fulspec(fk(:,:,j))));
   end
   
   return
   
%-------------------------------------------------------------------

function fkf = fulspec(fk);

    %     Assumes fk contains upper-half plane of spectral field, 
    %     and specifies lower half plane by conjugate 
    %     symmetry.  Input fk should be (1:2*kmax+1,1:kmax+1), kmax = 2^n-1.
    %     fkf is padded with zeros on (1,:) and (:,1), as expected
    %     by fftshift.  Grid resolution will be 2^(n+1) x 2^(n+1).  

    global kmax nx nkx nky

    fkf = zeros(nx,nx);
    fup = fk;
    fup(kmax:-1:1,1) = conj(fup(kmax+2:nkx,1));
    fdn = conj(fup(nkx:-1:1,nky:-1:2));
    fkf(2:nx,nky+1:nx) = fup;
    fkf(2:nx,2:nky) = fdn;
        
    return
    
%-------------------------------------------------------------------

function [hc] = plotstuff(frame,particles,ep,gam)
    
    persistent cvec
    global h zeta divuk xp yp nx t
    
    %if (frame==1)
    %   figure(10)
    %   axis square
    %   cvec = [min(h(:)) max(h(:))]
    %   disp('move and reshape figure 10 as desired, then press any key')
    %   pause
    %end
    
    % linear PV = gam*zeta-h
    
    figure(11);
    clf
    %plot(real(h(:,1)))
    %axis([1 size(h,1) -.1 .1])
    x = linspace(0,2*pi*(nx-1)/nx,nx);
    pcolor(x,x,ep*(gam*real(zeta)'-real(h)')), shading interp, colorbar, axis image
    title(strcat('PV at t=',num2str(t)),'fontsize',16)
    if (particles)
        hold
        %plot(xp,yp,'k.')
        plot(mod(xp,2*pi),mod(yp,2*pi),'k.')
    end
        
    %pcolor(real(zeta)'), shading interp, colorbar, axis tight manual
    %pcolor(k2g(divuk)'), shading interp, colorbar, axis tight manual
    %caxis(cvec)
    set(gca,'fontsize',16)

    drawnow
    hc=gca;

    return
