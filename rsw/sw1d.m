function [U,KE,PE,time] = sw1d(Ui,Ro,Bu,V0,Kv,nt,savestep)

%  [U,KE,PE,time] = sw1d(Ui,Ro,Bu,V0,Kv,nt,savestep)
%
%  Solve forced 1D Rotating Shallow Water equations 
%
%      u_t - v + Bu h_x + Ro uu_x = 0
%      v_t + u + Ro uv_x          = -Ro uV_x   
%      h_t + u_x + Ro (hu)_x      = 0
%
%  where H = 1+eps*h is fluid thickness, Ro = U/fL, Bu = g'H1/(fL)^2
%  and V_x(x) = - V0*sin(Kv*x) is an imposed periodic barotropic velocity
%  [associated with a free surface variation occuring on the veryasdkj
%  large external deformation scale (g H0)^(1/2)/f.]
%
%  Uses pseudospectral method and RK3 timestepping.  Inputs are
%  initial values of u, v, h in array Ui, with 
%      Ui(:,1) = u, Ui(:,2) = v, Ui(:,3) = h
%  Size of Ui is [nx,3], and nx should be a power of 2.
%  Additional inputs are Rossby number (Ro), number of timesteps
%  (nt) and timestep inteval at which to compute and save energy
%  and output fields.  Size of U is [nx,3,floor(nt/savestep)].

% Defaults
dttune= .1;

global KMAX KMAXBIG K NX NXBIG RO BU VX

NX = size(Ui,1); 
KMAX = NX/2 - 1;              % NX is a power of 2
KMAXBIG = 3*(KMAX+1)/2-1;     % Max k for padded spectral function
NXBIG = 2*(KMAXBIG+1);        % grid res for padded transform
K = (0:KMAX)';
RO = Ro;
BU = Bu;

% Constant coefficients for RK3 timestepping
c1 = 1/3; c2 = 5/9; c3 = 15/16; c4 = 153/128; c5 = 8/15;

% Get inital spectral fields
uk = g2s(Ui);

% Get forcing flow
if (V0~=0)
  x = linspace(0,2*pi,NXBIG+1); x=x(1:end-1);
  VX = -V0 * sin(Kv*x');
else
  VX = 0;
end

% for adaptive dt
cgw = sqrt(Bu + 1);  % assume k=1
dtfac = dttune*2*pi/KMAX;

frame = 0;
t=0;
for n = 1:nt

  % set dt
  vmax = max(cgw,sqrt(max(abs(Ui(:,1)))^2+max(abs(Ui(:,2)))^2));
  dt = dtfac/vmax;

  % Diagnostics and saves
  if ((mod(n,savestep)==0)||(n==1))
    frame = frame+1;
    U(:,:,frame) = s2g(uk);
    H = 1+Ro*U(:,3,frame);
    KE(frame) = sum(.5*H.*(U(:,1,frame).^2+U(:,2,frame).^2));
    PE(frame) = sum(.5/Ro^2 * H.^2);
    time(frame) = t;
    %dto(frame) = dt;
    %msg('KE =',KE(frame))
    %msg('dt=',dt)
  end
  
  % RK3 timestep
  rk  = dt*rhs(uk);
  uk1 = uk  + c1*rk;
  rk1 = dt*rhs(uk1) - c2*rk;
  uk2 = uk1 + c3*rk1;
  uk  = uk2 + c5*(dt*rhs(uk2) - c4*rk1);
  
  t = t+dt;
  
end

return 

%............................................................

function rk = rhs(uk)

global KMAX KMAXBIG K NXBIG RO BU VX

% KMAX = 2^n - 1
% NX = 2*(KMAX+1) = 2^(n+1) 
% KMAXBIG = 3*(KMAX+1)/2-1 = 3*2^(n-1)-1
% NXBIG = 2*(KMAXBIG+1) = 3*2^n
%
% Standard format for unshifted spectrum is K = [0:KMAX -KMAX-1:-1], 
% so starting from uk on K = [0:KMAX], first pad to 
% K = [0:KMAXBIG], then use conj symmetry to get it on 
% K = [0:KMAXBIG -KMAXBIG-1:-1], do ifft, nonlin terms, fft.
% Note that value at K=-KMAXBIG-1 is set to 0, and K=0 value
% (uk(1)) is not repeated 

wk = zeros(size(uk));
% We need products u*u, u*v_x and u*h, so first make wk s.t. 
wk(:,1) = uk(:,1);      % u
wk(:,2) = i*K.*uk(:,2);  % v_x
wk(:,3) = uk(:,3);      % h

wkbig = zeros(KMAXBIG+1,3);  
wkbig = [wk; zeros((KMAX+1)/2,3)];  % uk padded to K=[0:KMAXBIG]
wkbigf = zeros(NXBIG,3);
wkbigf(1:KMAXBIG+1,:) = wkbig;   
wkbigf(KMAXBIG+3:end,:) = conj(wkbig(end:-1:2,:)); % by conj symm
wgbig = real(ifft(wkbigf));  % *NXBIG for scaled result

% Now get nonlinear terms in x-space, then go back to k-space

wnl = zeros(NXBIG,3);
wnl(:,1) = wgbig(:,1).*wgbig(:,1);         % u*u
wnl(:,2) = wgbig(:,1).*(wgbig(:,2)+VX);    % u*(v_x+V_x)
wnl(:,3) = wgbig(:,1).*wgbig(:,3);         % u*h

% if we'd scaled ifft, all above products would be NXBIG^2 larger
% fft below should be divided by NXBIG, therefore result must be
% multiplied by NXBIG

wnlkbigf = NXBIG*fft(wnl);  
wnlk = wnlkbigf(1:KMAX+1,:);

% Finally compute RHS terms w/ <.> as transform 
% du/dt =  v -i*K*h -Ro*i*K*(1/2)*<u^2>
% dv/dt = -u        -Ro*<u*v_x>
% dh/dt = -i*K*u    -Ro*i*K*<h*u>

rk(:,1) =  uk(:,2) -BU*i*K.*uk(:,3) - .5*RO*i*K.*wnlk(:,1); 
rk(:,2) = -uk(:,1)               -    RO*wnlk(:,2);
rk(:,3) = -i*K.*uk(:,1)          -    RO*i*K.*wnlk(:,3);

%............................................................

function fg = s2g(fk)

% k-space to grid space

global KMAX NX

fkf = zeros(NX,3);
fkf(1:KMAX+1,:) = fk;   
fkf(KMAX+3:end,:) = conj(fk(end:-1:2,:)); % by conj symm
fg = NX*real(ifft(fkf));

%............................................................

function fk = g2s(fg)

% grid space to k-space

global KMAX KMAXBIG NX

fk = (1/NX)*fft(fg);   % on K=[0:KMAX -KMAX-1:-1]
fk = fk(1:KMAX+1,:);   % on K=[0:KMAX]
