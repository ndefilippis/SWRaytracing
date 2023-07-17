function [Sout,time,ke,pe,hmov] = swks(Sin,f,Cg,numsteps,savestep)

%  [Sout,time,ke,pe,hmov] = swks(Sin,f,Cg,numsteps,savestep)
%
%  Solves rotating shallow water equations (H = 1+h)
%
%  u_t =  v*(f + zeta) - B_x + nu*del^(a) u
%  v_t = -u*(f + zeta) - B_y + nu*del^(a) v
%  h_t = -[u*(h+1)]_x - [v*(h+1)]_y 
%  
%  w/ zeta = v_x-u_y and B = (u^2+v^2)/2 + Cg^2*h    
%
%  Inputs
%
%  Sin:       N x N x 3 array containing initial u, v, h, respectively
%  f:         Nondim Coriolis [ie inverse Rossby number f_0*L/U]
%  Cg:        Nondim GW speed  [sqrt(g*H_0)/U]
%  numsteps:  Total number of timesteps 
%  savestep:  Frequency, in timesteps, to save output    
%
%  Outputs
%
%  Sout:      Arranged as Sin, but with 4th dimension for time
%  time:      Times at which output is saved
%  ke:        Time series of KE
%  pe:        Time series of PE
%  hmov:      Movie of h field (if hmov included in output list)
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
dttune = .5;
nutune = .2;

% Check for outputs requested
makemov=false;
if (nargout>4), makemov=true;  end
   
% Set global params for use in rhs functions
global nx kgfac gkfac1 gkfac2 ikx_ iky_ kmax nkx nky u v h zeta divuk

% Get and check dimensions
[nx,ny] = size(Sin(:,:,1));
if (nx~=ny), error('must have nx = ny'); end
if (mod(log2(nx),1)~=0), error('must have nx = 2^n, n integer'); end

% Dimensions of spectral fields -- stored on upper half-plane in kx,ky
% (lower half-plane given by conjugate symmetry).
kmax = nx/2 - 1;  % 2^n - 1
nkx = 2*kmax+2;   % -1-kmax:kmax, nkx = nx
nky = kmax+1;     % 0:kmax
dx = 2*pi/nx;     % L_domain = 2*pi

% Set up arrays of wavenumbers for computing derivatives
[ikx_,iky_] = ndgrid([0:kmax -1-kmax:-1],0:kmax);  % fft order
K_   = sqrt(ikx_.^2 + iky_.^2);
ikx_ = sqrt(-1)*ikx_; 
iky_ = sqrt(-1)*iky_;

% Initialize fields for spectral products with Orszag dealiasing
kcut      = sqrt(8./9.)*(kmax+1);
damask    = ones(size(ikx_));
damask(K_>kcut)    = 0.;
damask(kmax+2:end,1)   = 0.;  % kx<0, ky=0 given by conj sym

eipik  = exp(pi*(ikx_+iky_)/nx);
kgfac  = (1 + sqrt(-1)*eipik).*damask; 
gkfac1 = (1 - sqrt(-1)*conj(eipik))/4;
gkfac2 = (1 + sqrt(-1)*conj(eipik))/4;

Ka_ = K_.^a;  % For hyperviscous filter

clear K_ eipik damask

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

% Set counters and temp fields
frame = 0;  t = 0;  n = 0; 
Rk = zeros(size(Sk)); Rkm1 = Rk; Rkm2 = Rk;

keepgoing = true;
while keepgoing
        
    % Save n-1 and n-2 rhs and get next one.  
    % getrhs sets u, v, h, too, but they have imaginary parts
    % holding staggered grid fields; use real() for computations etc 
    Rkm2 = Rkm1;
    Rkm1 = Rk;
    Rk = getrhs(Sk,f,Cg);  
    
    if (n==0) Rkm1 = Rk; Rkm2 = Rk; end

    Umax = max([max(max(abs(real(u)))) max(max(abs(real(v)))) Cg]);

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
        ke(frame) = .5*sum(real(1+h(:)).*(real(u(:)).^2+real(v(:)).^2));
        pe(frame) = .5*Cg^2*sum(real(1+h(:)).^2);
        %hbar(frame) = sum(real(h(:)))/nx^2;
        Sout(:,:,1,frame) = real(u);    
        Sout(:,:,2,frame) = real(v);    
        Sout(:,:,3,frame) = real(h);  
        time(frame) = t;
        disp(strcat('Wrote frame >',num2str(frame),' out of >',num2str(nframes)))
        disp(strcat('max(|u|) = ',num2str(Umax),', dt = ',num2str(dt),', nu = ',num2str(nu)))
        if (makemov)
            hc = plotstuff(frame);
            hmov(frame) = getframe(hc);
        end
    end
    
    % Timestep and diffuse
    Sk = repmat(filterdn,[1 1 3]).*(repmat(filterup,[1 1 3]).*Sk + dt*(a1*Rk + a2*Rkm1 + a3*Rkm2));

    n = n+1;
    t = t+dt;  % clock
    
    if (n==numsteps), disp('End reached'), keepgoing=false; end

end

return

%-------------------------------------------------------------------
% Internal functions
%-------------------------------------------------------------------

function Rk = getrhs(Sk,f,Cg)
   
    global ikx_ iky_ u v h zeta divuk

    u = k2gp(Sk(:,:,1));
    v = k2gp(Sk(:,:,2));
    h = k2gp(Sk(:,:,3));

    zeta = k2gp(ikx_.*Sk(:,:,2) - iky_.*Sk(:,:,1));     % vorticity
    divuk = ikx_.*Sk(:,:,1) + iky_.*Sk(:,:,2);
    Bk = gp2k(gprod(u,u)+gprod(v,v)) + Cg^2*Sk(:,:,3);  % Bernouli in k-spc
    
    Rk(:,:,1) =  gp2k(gprod(v,zeta)) + f*Sk(:,:,2) - ikx_.*Bk;
    Rk(:,:,2) = -gp2k(gprod(u,zeta)) - f*Sk(:,:,1) - iky_.*Bk;
    Rk(:,:,3) = -ikx_.*gp2k(gprod(u,h)) - iky_.*gp2k(gprod(v,h)) -divuk;
    
    return
    
%-------------------------------------------------------------------
        
function fg = k2gp(fk)

   % Transform to grid space, packing fields shifted by dx/2 into
   % imaginary part
    
   global nx kgfac

   fkt = fulspec(kgfac.*fk);
   fg  = nx^2*ifft2(fkt);

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
  
    Wk  = fft2(prodg)/nx^2;
    
    % Extract spectral products on grid and shifted grid, and average.
    Wk_up = Wk(:,1:kmax+1); % uhp, ordered as ikx_ etc...
    Wk_dn = zeros(size(Wk_up));
    Wk_dn(1,1) = conj(Wk(1,1));
    Wk_dn(2:kmax+1,1) = conj(Wk(end:-1:kmax+3,1));
    Wk_dn(kmax+3:end,1) = conj(Wk(kmax+1:-1:2,1));
    Wk_dn(1,2:kmax+1) = conj(Wk(1,end:-1:kmax+3));
    Wk_dn(2:kmax+1,2:kmax+1) = conj(Wk(end:-1:kmax+3,end:-1:kmax+3));
    Wk_dn(kmax+3:end,2:kmax+1) = conj(Wk(kmax+1:-1:2,end:-1:kmax+3));
    
    prodk = gkfac1.*Wk_up + gkfac2.*Wk_dn;
    
    return

%-------------------------------------------------------------------
    
function fk = g2k(fg)
    
    % Just for transforming gridded input field
        
    global nx kmax
    
    fk = zeros(nx,kmax+1,3);
    
    for j=1:size(fg,3)
        fkt = fft2(fg(:,:,j))/nx^2;
        fk(:,:,j) = fkt(:,1:kmax+1);  % keep only ky>=0 part      
    end
    
    return
    
%-------------------------------------------------------------------
        
function fg = k2g(fk)

   % Transform to grid space - just for output.
    
   global nx 

   fg = zeros(nx,nx,3);
   
   for j=1:size(fk,3)
       fg(:,:,j)  = nx^2*ifft2(fulspec(fk(:,:,j)));
   end
   
   return
   
%-------------------------------------------------------------------

function fkf = fulspec(fk);

    %     Assumes fk contains upper-half plane of spectral field, 
    %     and specifies lower half plane by conjugate 
    %     symmetry. Input has size [nx nx/2], arranged as 
    %     [0:kmax -1-kmax:-1,0:kmax] .. recall nx = 2*kmax+2
    %     kx=-1-kmax and ky=-1-kmax are padded with 0's
        
    global kmax nx nkx nky

    fkf = zeros(nx,nx);
    fkf(1:end,1:kmax+1) = fk;
    fkf(kmax+3:end,1) = conj(fkf(kmax+1:-1:2,1));
    fkf(1,kmax+3:end) = conj(fkf(1,kmax+1:-1:2));
    fkf(2:kmax+1,kmax+3:end) = conj(fkf(end:-1:kmax+3,kmax+1:-1:2));
    fkf(kmax+3:end,kmax+3:end) = conj(fkf(kmax+1:-1:2,kmax+1:-1:2)); 
        
    return
    
%-------------------------------------------------------------------

function [hc] = plotstuff(frame)
    
    persistent cvec
    global h zeta divuk
    
    %if (frame==1)
    %   figure(10)
    %   axis square
    %   cvec = [min(h(:)) max(h(:))]
    %   disp('move and reshape figure 10 as desired, then press any key')
    %   pause
    %end
    
    figure(10);
    clf
    pcolor(1+real(h)')
    %pcolor(real(zeta)')
    %pcolor(k2g(divuk)')
    shading interp
    colormap(jet)
    %caxis(cvec)
    colorbar
    axis tight manual
    title('h')
    set(gca,'fontsize',16)
    drawnow
    hc=gca;

    return
