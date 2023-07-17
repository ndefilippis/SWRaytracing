function [Sout,time,ke,pe,Xp,hmov] = swkUqx(Sin,Psi,f,Cg,numsteps,savestep,Xpin)

%  [Sout,time,ke,pe,Xp,hmov] = swk(Sin,Psi,f,Cg,numsteps,savestep,Xpin)
%
%  Solves rotating shallow water equations (H = 1+h)
%
%  u_t = -U*u_x - V*u_y - u*U_x - v*U_y + f*v - Cg^2*h_x + nu*del^(a) u
%  v_t = -U*v_x - V*v_y - u*V_x - v*V_y - f*u - Cg^2*h_y + nu*del^(a) v
%  h_t = -(u*h)_x -(v*h)_y -(u_x + v_y) 
%
%  where U = - Psi_y, V = Psi_x 
%
%  Inputs
%
%  Sin:       N x N x 3 array containing initial u, v, h, respectively
%  f:         Nondim Coriolis [ie inverse Rossby number f_0*L/U]
%  Cg:        Nondim GW speed  [sqrt(g*H_0)/U]
%  numsteps:  Total number of timesteps 
%  savestep:  Frequency, in timesteps, to save output    
%  Xpin:      Structure Xpin.x, Xpin.y containing intial particle 
%             positions (optional)
%
%  Outputs
%
%  Sout:      Arranged as Sin, but with 4th dimension for time
%  time:      Times at which output is saved
%  ke:        Time series of KE
%  pe:        Time series of PE
%  hmov:      Movie of h field (if hmov included in output list)
%  Xp:        Structure with coordinates Xp(j).x, Xp(j).y of particles,
%             where j is timestep
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
nutune = 1;
dttune = .1;   % Courant number

pvdamprate = 1.;

% Set global params for use in rhs functions
global nx ikx_ iky_ kmax nkx nky u v h zeta divuk U V
global kgfac gkfac1 gkfac2 damask 
global t


if nargin > 6, 
    particles=1; 
    Xp = Xpin;
    np = length(Xpin.x);
else 
    particles=0;
    Xp=0;
end

% Check for outputs requested
makemov=false;
if (nargout>5), makemov=true;  end

% Maximum phase speed.  Note omega = +/-sqrt(Cg^2*K^2+f^2) 
Cmax = sqrt(Cg^2+f^2);

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
ikx_ = 1i*ikx_; 
iky_ = 1i*iky_;

% Initialize fields for spectral products with Orszag dealiasing
kcut      = sqrt(8./9.)*(kmax+1);
damask    = ones(size(ikx_));
damask(K_>kcut)    = 0.;
damask(1:kmax,1)   = 0.;  % changed from 1:kmax+1 to keep 0,0 value.

eipik  = exp(pi*(ikx_+iky_)/nx);
kgfac  = 1 + 1i*fulspec(eipik);
gkfac1 = (1 - 1i*conj(eipik))/4;
gkfac2 = (1 + 1i*conj(eipik))/4;

% Trapezoidal hyperdiffusion operators
nudt = nutune*2*pi/(nx*kmax^a);  % divide by dt to get effective nu
fR = (1+nudt/2*K_.^a).^(-1);
fU = (1-nudt/2*K_.^a).*fR;
filterU = ones(nkx,nky,3);
filterR = ones(nkx,nky,3);
filterU(:,:,1:2) = repmat(fU,[1 1 2]);       
filterR(:,:,1:2) = repmat(fR,[1 1 2]);
clear K_ f1 f2

% Get initial spectral PV
Sk = g2k(Sin);

Psik = g2k(Psi);
U = -k2g(iky_.*Psik);
V = k2g(ikx_.*Psik);

% Set params for AB3+trapezoidal diffusion (Durran 3.81)
a1 = 23/12;   a2 = -16/12;  a3 = 5/12;  

% Preallocate arrays to hold saved output
nframes = floor(numsteps/savestep);
Sout = zeros(nx,ny,3,nframes);  % output on grid
time = zeros(1,nframes); 
ke = zeros(1,nframes); 
pe = zeros(1,nframes); 

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
    Rk = getrhs(Sk,f,Cg);
    
    if (n==0) Rkm1 = Rk; Rkm2 = Rk; end

    Umax = max([max(abs(real(u(:)))) max(abs(real(v(:)))) Cmax]);

    % Exit if blowing up, and save last field.
    if (Umax>1e6)  
        disp(strcat('Blow up! Umax= ',num2str(Umax),', t= ',num2str(t)))
        Sout(:,:,:,frame+1) = k2g(Sk); 
        keepgoing = false;
    end
    
    % Adapt dt and nu
    dt = dttune*dx/Umax;   % Courant condition 
    nu = nudt/dt;
    
    % Save output at frequency savestep
    if (mod(n,savestep)==0||n==0)  
        frame = frame+1;  
        ke(frame) = .5*sum(real(1+h(:)).*(real(u(:)).^2+real(v(:)).^2))/nx^2;
        pe(frame) = .5*Cg^2*sum(real(h(:)).^2)/nx^2;
        %hbar(frame) = sum(real(h(:)))/nx^2;
        Sout(:,:,1,frame) = real(u);    
        Sout(:,:,2,frame) = real(v);    
        Sout(:,:,3,frame) = real(h);  
        time(frame) = t;
        disp(strcat('Wrote frame >',num2str(frame),' out of >',num2str(nframes)))
        disp(strcat('max(|u|) = ',num2str(Umax),', dt = ',num2str(dt),', nu = ',num2str(nu)))
        if (makemov)
            if (particles)
                xp=[Xp(n+1).x]; yp=[Xp(n+1).y];
                hc = plotstuff(frame,xp,yp);
            else
                hc = plotstuff(frame);
            end            
            hmov(frame) = getframe(hc);
        end
        %        if writeoutput
        %            save(outputfile,'Sout','ke','pe','time','Xp')
        %        end
        save
    end
    
    % Timestep and diffuse
    Sk = filterU.*Sk + dt*filterR.*(a1*Rk + a2*Rkm1 + a3*Rkm2);
    
    Sk = PV_damping(Sk);

    if (particles) % advect particles
        Xp(n+2) = advect_particles(Xp(n+1),real(u),real(v),dx,dx,dt);
    end
    
    n = n+1;
    t = t+dt;  % clock
    
    if (n==numsteps), disp('End reached'), keepgoing=false; end

end

return

%-------------------------------------------------------------------
% Internal functions
%-------------------------------------------------------------------

function Rk = getrhs(Sk,f,Cg)
   
%  u_t = -U*u_x - V*u_y - u*U_x - v*U_y + f*v - Cg^2*h_x + nu*del^(a) u
%      = -(U*u)_x - (V*u)_y + U*divu + f*v - Cg^2*h_x + nu*del^(a) u
%  v_t = -U*v_x - V*v_y - u*V_x - v*V_y - f*u - Cg^2*h_y + nu*del^(a) v
%      = -(U*v)_x - (V*v)_y + V*divu - f*u - Cg^2*h_y + nu*del^(a) v
%  h_t = -(u*h)_x -(v*h)_y -(u_x + v_y) 

    global ikx_ iky_ u v U V h zeta divuk

    u = k2gp(Sk(:,:,1));
    v = k2gp(Sk(:,:,2));
    h = k2gp(Sk(:,:,3));

    zeta = k2gp(ikx_.*Sk(:,:,2) - iky_.*Sk(:,:,1));     % vorticity
    divuk = ikx_.*Sk(:,:,1) + iky_.*Sk(:,:,2);
    divu = k2gp(divuk);
    
    Rk(:,:,1) = - 2*ikx_.*gp2k(gprod(U,u)) ...
                -   iky_.*gp2k( gprod(V,u) + gprod(v,U )) ...
                +         gp2k(gprod(U,divu)) ...
                +  f*Sk(:,:,2) - Cg^2*ikx_.*Sk(:,:,3);
    Rk(:,:,2) = -  ikx_.*gp2k(gprod(U,v)+ gprod(u,V)) ...
                -2*iky_.*gp2k(gprod(V,v)) ...
                +        gp2k(gprod(V,divu)) ...
                - f*Sk(:,:,1) - Cg^2*iky_.*Sk(:,:,3);
    Rk(:,:,3) = -ikx_.*gp2k(gprod(U,h)) -iky_.*gp2k(gprod(V,h))- divuk;
    
    return

    
%-------------------------------------------------------------------
        
function Sk_new = PV_damping(Sk,f,Cg,pvdamprate)
   global ikx_ iky_ 

   zeta = k2gp(ikx_.*Sk(:,:,2) - iky_.*Sk(:,:,1)); 
   h = k2gp(Sk(:,:,3));
   pv_res = zeta - h*f;
   
   psi_res = gp2k(pv_res)./(Cg^2/f*ikx_.^2 + Cg^2/f*iky_.^2 - f);
   u_res = - iky_.*psi_res*Cg^2/f;
   v_res = ikx_.*psi_res*Cg^2/f;
%    psi_res = gp2k(pv_res)./(ikx_.^2 + iky_.^2);
%    psi_res(128,1) = 0;
%    u_res = - iky_.*psi_res;
%    v_res = ikx_.*psi_res;

   Sk_new = Sk;
   Sk_new(:,:,1) = Sk_new(:,:,1) - u_res*pvdamprate;
   Sk_new(:,:,2) = Sk_new(:,:,2) - v_res*pvdamprate;
   Sk_new(:,:,3) = Sk_new(:,:,3) - psi_res*pvdamprate;
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
        
    prodg = real(f).*real(g) + 1i*imag(f).*imag(g);

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

function [hc] = plotstuff(frame,xp,yp)
    
    persistent cvec x
    global h zeta divuk nx t 
    
    if (frame==1)
        x = linspace(0,2*pi*(nx-1)/nx,nx)-pi;
        figure(11)
        clf
        axis square
        disp('move and reshape figure 11 as desired, then press any key')
        pause
    end

    pcolor(x,x,real(h)'), shading interp, colorbar, axis image
    clim = max([abs(min(real(h(:)))) abs(max(real(h(:))))]);
    caxis([-clim clim])

    %plot(x,real(h(:,1))), grid
    %axis([-pi pi -.2 .2])
    
    %subplot(1,2,1)
    %pcolor(x,x,real(zeta-h)'), shading interp, colorbar, axis image
    %title(strcat('PV at t=',num2str(t)),'fontsize',14)
    %if (nargin>1)
    %    hold
    %    plot(xp,yp,'k.','MarkerSize',10)
    %end
    %set(gca,'fontsize',14)
    %subplot(1,2,2)
    %pcolor(x,x,k2g(divuk)'), shading interp, colorbar, axis image
    %title('Divergence','fontsize',14) 
    %set(gca,'fontsize',14)
    
    %pcolor(real(zeta)'), shading interp, colorbar, axis tight manual
    %pcolor(k2g(divuk)'), shading interp, colorbar, axis tight manual
    %caxis(cvec)
   

    drawnow
    hc=gca;

    return
