function [A,time] = ybj1d(A0,Bu,V0,Kv,dt,nt,savestep)

%  [A,t] = ybj1d(A0,Bu,V0,Kv,dt,nt,savestep)
%
%  Solve forced 1D YBJ amplitude equation
%
%      A_T + (i/2)(V_x A - Bu A_xx) = 0
%
%  where V_x(x) = -V0*sin(Kv*x) is an imposed periodic barotropic
%  vorticity. 
%  
%  Uses pseudospectral method and RK3 timestepping.  Inputs are
%  initial complex A0 = A(x,0), with size nx, where nx should be a
%  power of 2.  Additional inputs are Burger number Bu, number of
%  timesteps (nt) and timestep at which to save output fields
%  (savestep).  Size of a is [nx,floor(nt/savestep)].


global KMAX KMAXBIG K NX NXBIG BU VX

NX = length(A0); 
KMAX = NX/2 - 1;              % NX is a power of 2
KMAXBIG = 3*(KMAX+1)/2-1;     % Max k for padded spectral function
NXBIG = 2*(KMAXBIG+1);        % grid res for padded transform
K=[0:KMAX -KMAX-1:-1]';
BU = Bu;

% Constant coefficients for RK3 timestepping
c1 = 1/3; c2 = 5/9; c3 = 15/16; c4 = 153/128; c5 = 8/15;

% Get inital spectral field ..  A is complex!
Ak = (1/NX)*fft(A0);  % on K=[0:KMAX -KMAX-1:-1]

% Get forcing field on big K grid
x = linspace(0,2*pi,NXBIG+1); x=x(1:end-1);
VX = -V0 * sin(Kv*x');

frame = 0;
t=0;
for n = 1:nt

  % Diagnostics and saves
  if ((mod(n,savestep)==0)||(n==1))
    frame = frame+1;
    A(:,frame) = NX*ifft(Ak);
    time(frame) = t;
  end
  
  % RK3 timestep
  rk  = dt*rhs(Ak);
  Ak1 = Ak  + c1*rk;
  rk1 = dt*rhs(Ak1) - c2*rk;
  Ak2 = Ak1 + c3*rk1;
  Ak  = Ak2 + c5*(dt*rhs(Ak2) - c4*Ak1);
  
  t = t+dt;
  
end

return 

%............................................................

function rk = rhs(Ak)

global KMAX KMAXBIG K NXBIG BU VX

% KMAX = 2^n - 1
% NX = 2*(KMAX+1) = 2^(n+1) 
% KMAXBIG = 3*(KMAX+1)/2-1 = 3*2^(n-1)-1
% NXBIG = 2*(KMAXBIG+1) = 3*2^n
%
% Standard format for unshifted spectrum is 
% K = [0:KMAX -KMAX-1:-1] = [1:KMAX+1 KMAX+2:NX],  
% so starting from Ak on K, first pad to 
% Kbig = [0:KMAXBIG -KMAXBIG-1:-1] do ifft, nonlin terms, fft.
% Note that value at K=-KMAXBIG-1 is set to 0, and K=0 value
% (Ak(1)) is not repeated 

Akbig = zeros(NXBIG,1);  
Akbig(1:KMAX+1) = Ak(1:KMAX+1);
Akbig(end-KMAX:end) = Ak(end-KMAX:end);
Abig = ifft(Akbig);  % *NXBIG for scaled result

% Now get nonlinear terms in x-space, then go back to k-space
AVXbig = Abig.*VX;   

% if we'd scaled ifft, all above products would be NXBIG^2 larger
% fft below should be divided by NXBIG, therefore result must be
% multiplied by NXBIG

AVXkbig = NXBIG*fft(AVXbig);
AVXk = zeros(size(Ak));
AVXk(1:KMAX+1) = AVXkbig(1:KMAX+1);
AVXk(end-KMAX:end) = AVXkbig(end-KMAX:end);


% Finally compute RHS terms w/ <.> as transform 
% dA/dt = (i/2)(VX*A  - Bu*A_xx)

rk = -(sqrt(-1)/2)*(AVXk + BU*K.^2.*Ak);

