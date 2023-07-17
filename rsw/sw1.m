function [Uout,ke,pe,time,Xp,t,ue] = sw1(Ui,x,f,Cg,numsteps,savestep,Xpin)

%  [Uout,ke,pe,time,Xp,t,ue] = sw1(Ui,x,f,Cg,numsteps,savestep,Xpin)
%
%  Solve 1D Rotating Shallow Water equations 
%
%      u_t - f v + Cg^2 h_x + .5(uu)_x = 0
%      v_t + f u            + uv_x     = 0  
%      h_t + u_x            + (hu)_x   = 0
%
%  where H = H0*(1+h) is fluid thickness and Cg = sqrt(g*H_0)
%
%  Inputs are initial values of u, v, h in array Ui, with 
%      Ui(:,1) = u
%      Ui(:,2) = v
%      Ui(:,3) = h
%  Size of Ui is [NX,3], and NX should be a power of 2.
%  Additional inputs are grid for the model x(:) (length NX),
%  number of timesteps (numsteps) and timestep inteval at which to compute and save energy
%  and output fields (savestep).  Optional Xpin(:) (length NP) are initial
%  positions for a set of NP Lagrangian particles.  
%  Size of output Uout(:,:,:) is [NX,3,floor(nt/savestep)].  Also gives
%  energies ke(:) and pe(:), and time(:) at which fields and
%  diagnostics are saved. If Xpin specified, output Xp(:,:) has
%  size (NP,numsteps).
%  Uses pseudospectral method, and AB3 timestepping with trapezoidal
%  diffusion (Durran 3.81).  Particle advection uses RK4

% Defaults
a = 8;         % hyperviscosity order, nu del^a u 
nutune = .01;
dttune = .3;

% For particle advection
if nargin>5
    np = length(Xpin);
    Xp = zeros(np,numsteps);
    Xp(:,1) = Xpin;  
    particles = 1;
else
    particles = 0;
end

global KMAX KMAXBIG K NX NXBIG CG2 F 

NX = size(Ui,1); 
KMAX = NX/2 - 1;              % NX is a power of 2
KMAXBIG = 3*(KMAX+1)/2-1;     % Max k for padded spectral function
NXBIG = 2*(KMAXBIG+1);        % grid res for padded transform
K = (0:KMAX)';
CG2 = Cg^2;
F = f;

dx = 2*pi/NX;                 % grid size
Cmax = sqrt(Cg^2 + f^2);      % IGW speed - assume k=1

% Set params for AB3+trapezoidal diffusion (Durran 3.81)
a1 = 23/12;   a2 = -16/12;  a3 = 5/12;  
Ka = K.^a;                    % For hyperviscous filter

Uk = g2s(Ui);                 % Get initial spectral fields
ue(1) = Uk(1,1);

% Preallocate arrays to hold saved output
nframes = floor(numsteps/savestep);
Uout = zeros(NX,3,nframes);  % output on grid
time = zeros(1,nframes); 
ke = zeros(1,nframes); 
pe = zeros(1,nframes); 
filterup = ones(KMAX+1,3);
filterdn = filterup;
ue = zeros(1,numsteps); % store mean u
% Set counters and temp fields
frame = 0;  t = 0;  n = 0; 
Rk = 0; Rkm1 = 0; Rkm2 = 0;

keepgoing = true;
while keepgoing

    n = n+1;
    
    U = s2g(Uk);
    
    % Diagnostics and saves
    if ((mod(n,savestep)==0)||(n==1))
        frame = frame+1;
        Uout(:,:,frame) = U;
        H = 1+U(:,3);
        ke(frame) = sum(.5*H.*(U(:,1).^2+U(:,2).^2));
        pe(frame) = sum(.5 * CG2 * H.^2);
        time(frame) = t(n);
        disp(strcat('Wrote frame >',num2str(frame),' out of >',num2str(nframes)))
        %disp(strcat('max(|u|) = ',num2str(vmax),', dt = ',num2str(dt),', nu = ',num2str(nu)))
        %if (makemov)
        %    hc = plotstuff(frame,particles);
        %    hmov(frame) = getframe(hc);
        %end
    end
    
    % Save n-1 and n-2 rhs and get next one.  
    % getrhs sets u, v, h, too, but they have imaginary parts
    % holding staggered grid fields; use real() for computations etc 
    Rkm2 = Rkm1;
    Rkm1 = Rk;
    Rk = rhs(Uk);  
    if (n==1) Rkm1 = Rk; Rkm2 = Rk; end

    % Check max vel for dt and 
    vmax = max([max(abs(real(U(:,1)))) max(abs(real(U(:,2)))) Cmax]);

    if (vmax>1e6)  % Exit if blowing up
        disp(strcat('Blow up! vmax= ',num2str(vmax),', t= ',num2str(t)))
        keepgoing = false;
    end

    % Adapt dt and nu
    dt = dttune*dx/vmax;   % Courant condition 
    nu = nutune*dx^a/dt;   % Optimal value of hyperviscous coef

    % Hyperviscous filters 
    filterup(:,1:2) = repmat(1 - (dt/2)*nu*Ka,[1 2]);       
    filterdn(:,1:2) = repmat(1./(1 + (dt/2)*nu*Ka),[1 2]);
    
    % Timestep and diffuse
    Uk = filterdn.*(filterup.*Uk + dt*(a1*Rk + a2*Rkm1 + a3*Rkm2));

    if (particles) % advect particles
        Xp(:,n+1) = advect1d(Xp(:,n),U(:,1),x,dt);
    end
    ue(n+1) = Uk(1,1);
    
    t(n+1) = t(n)+dt;  % clock
    
    if (n==numsteps), disp('End reached'), keepgoing=false; end

end

return 

%............................................................

function rk = rhs(uk)

global KMAX KMAXBIG K NXBIG CG2 F 

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

%UBIG = NXBIG*wgbig;  % For use in setting timestep etc...

% Now get nonlinear terms in x-space, then go back to k-space

wnl = zeros(NXBIG,3);
wnl(:,1) = wgbig(:,1).*wgbig(:,1);         % u*u
wnl(:,2) = wgbig(:,1).*wgbig(:,2);         % u*v_x
wnl(:,3) = wgbig(:,1).*wgbig(:,3);         % u*h

% if we'd scaled ifft, all above products would be NXBIG^2 larger
% fft below should be divided by NXBIG, therefore result must be
% multiplied by NXBIG

wnlkbigf = NXBIG*fft(wnl);  
wnlk = wnlkbigf(1:KMAX+1,:);

% Finally compute RHS terms w/ <.> as transform 
% du/dt =  f*v -i*Cg^2*K*h -i*K*(1/2)*<u^2>
% dv/dt = -f*u        -<u*v_x>
% dh/dt = -i*K*u    -i*K*<h*u>

rk(:,1) =  F*uk(:,2)    -CG2*i*K.*uk(:,3) -.5*i*K.*wnlk(:,1); 
rk(:,2) = -F*uk(:,1)                      -wnlk(:,2);
rk(:,3) = -i*K.*uk(:,1)                   -i*K.*wnlk(:,3);

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
