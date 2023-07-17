function [Diag_out,Params_out,F_out] = swp(F_in,varargin)

%  [Diag_out,Params_out,F_out] = swp(F_in,'parameter1',value1,'parameter2',value2,...)
%
%  Rotating Shallow Water model.  Solves equations:
%  
%  U_t + (f + zeta)*(u*jhat-v*ihat) = - grad(B) + Nu*(dx^2/dt)*del^2(U) - Dr*U
%  H_t + del.(UH) = Q
%
%  where 
%  U = u*ihat + v*jhat            is the vector velocity field
%  f = Roi + Be*y                 is the Coriolis parameter
%  zeta = dv/dx-du/dy             is the vorticity
%  B = (u^2+v^2)/2+Cg^2*h         is the Bernouli function
%  H = h - h_b                    is the fluid depth at each point
%  h_b                            is the bottom relief
%
%  Input fields -- initial conditions (required)
%    F_in:  structure containing initial u, v, h, optional topography h_b.
%    Set any field to 0 to initialize to 0.  
%    Structure must have fields with following sizes:
%    F_in.u:  (nx,ny-1) or 0   
%    F_in.v:  (nx-1,ny) or 0   
%    F_in.h:  (nx-1,ny-1) (required!)
%    F_in.hb: (nx-1,ny-1) or 0 
%    Additionally, you may specify F_in.frame (frame of input field) and 
%    F_in.time (model time of input fields) -- useful if they are restart fields
%  
%  Input parameters: 
%    Inputs are entered in the form "...,'parameter',value,..." in argument
%    list.  Default values are noted in brackets [].  
%
%    Physical parameters:
%      Roi:           inverse Rossby number {f_0*L/U_0} [0]
%      Beta:          beta parameter {beta_0*L^2/U_0}  [0]
%      Cg:            nondimen GW speed  (sqrt(g*H_0)/U_0) [0]
%      Drag:          nondimensional linear drag coef {r*L/U_0} [0]
%      Nu:            tuning factor for adaptive viscosity  [0]
%      Hdot:          nondimensional mass forcing {Hdotd*L/(H_0*U_0)} [0]
%      periodx:       flag to set for periodic BC in x [false]
%      periody:       flag to set for periodic BC in y [false]
%      geovel:        flag to set u, v to geostrophic vals [false]
%
%    Model control parameters:
%      nt:            number of timesteps to run the model  [100]
%      dttune:        tuning factor for adaptive timestep [0.2]
%      savestep:      timestep interval at which model state is saved [100]
%      diagstep:      timestep interval at which diagnostics done [savestep]
%      idstring:      string to append to ouput filenames  [''] 
%      datadir:       directory to which data will be written [./]
%      writetofiles:  flag to write fields to files  [false]
%                     (writes files u,v,h,zeta,q,time with idstring appended)
%                     If false, all frames written to F_out (access each
%                     field via, e.g., F_out(frame).v(:,:);).  If true,
%                     F_out contains only last frame for use as restart
%                     file.
%
%  Outputs:
%    Diag_out:  structure with fields
%      ke:   kinetic energy time series
%      ape:  available potential energy time series
%      htot: sum of fluid thickness time series
%      time: time at which diagnostics were written
%
%    F_out:  structure like F_in containing same fields, and additionally:
%      time: time at which fields were written
%
%  Numerics:
%  
%    Model uses RK3 and a C-grid with centered differencing. On
%    C-grid, boundaries run through vorticity points (i,j), h and B
%    points are at (i+1/2,j+1/2), u points at (i,j+1/2), v points at
%    (i+1/2,j).  Thus no normal flow conditions are applied explcitly,
%    since E-W boundaries run through u points, and similarly for v
%    points.  BC are either free-slip or periodic.  Time step is set
%    adaptively, multiplied by nondimensional prefactor dttune.
%    Calculated values of nu and dt are output to screen.

% Have to decide in avg and dif whether to assume point just beyond
% boundary is same to last interior point, or zero.  For vorticity, want
% zero.  For h, want h_(N+1) = h_N

%  For nonperiodic BC,
%  have problem with right and upper boundaries.  Waves emanate from those
%  Bs.  Probably with shift??  and zeros installed at those edges.

%% 3/18/05:  periodic BCs work in both directions!

%% UGH what did you do to lateral grid sizes? input struct is
% converted to array with same nx, ny for each field. ....

% Check for input parameters and set to defaults if not set by user
params = struct(varargin{:});
pf = fieldnames(params);
if (isempty(strmatch('nt',pf))),           params.nt = 500; end
if (isempty(strmatch('savestep',pf))),     params.savestep = 100; end 
if (isempty(strmatch('diagstep',pf))),     params.diagstep = params.savestep; end
if (isempty(strmatch('writetofiles',pf))), params.writetofiles = false; end
if (isempty(strmatch('idstring',pf))),     params.idstring = ''; end
if (isempty(strmatch('datadir',pf))),      params.datadir = ''; end
if (isempty(strmatch('Roi',pf))),          params.Roi = 0; end
if (isempty(strmatch('Beta',pf))),         params.Beta = 0; end
if (isempty(strmatch('Cg',pf))),           params.Cg = 0; end
if (isempty(strmatch('Drag',pf))),         params.Drag = 0; end
if (isempty(strmatch('Hdot',pf))),         params.Hdot = 0; end
if (isempty(strmatch('Nu',pf))),           params.Nu = 0; end
if (isempty(strmatch('periodx',pf))),      params.periodx = false; end
if (isempty(strmatch('periody',pf))),      params.periody = false; end
if (isempty(strmatch('geovel',pf))),       params.geovel = false; end
if (isempty(strmatch('dttune',pf))),       params.dttune = .2; end

% Check for input errors
if (~isstruct(F_in)), error('F_in not a structure'), end
Fnames = fieldnames(F_in);
if (~isempty(strmatch('frame',Fnames)))
    frame = F_in.frame;
else
    frame = 0;
end
if (~isempty(strmatch('time',Fnames))),        
    time = F_in.time;
else
    time = 0;
end

if (params.Beta~=0)&(params.periody), error('Cant have beta-plane with periodic-in-y'), end
if (params.Roi==0)&(params.geovel), error('Need f~=0 for geostrophic velocities'), end

% Set file names of files to be written
if (params.writetofiles)
    hname = strcat(params.datadir,'h',params.idstring);
    qname = strcat(params.datadir,'q',params.idstring);
    zetaname = strcat(params.datadir,'zeta',params.idstring);
    uname = strcat(params.datadir,'u',params.idstring);
    vname = strcat(params.datadir,'v',params.idstring);
    fname = strcat(params.datadir,'F_out',params.idstring);
end

% Height field
h = F_in.h;

% Extract dimensions
nx = size(h,1)+1;
ny = size(h,2)+1;

% Topography
if (~isempty(strmatch('hb',Fnames)))
    hb = F_in.hb;
    if (size(nb,1)~=nx-1&size(nb,2)~=ny), error('Size of hb wrong'), end
else
    hb = 0;
end

% Length scale
L = 1;

% Differentials
DX = L/nx;
DY = L/ny;
dr = 2*DX*DY/(DX+DY);

% Parameters
Cg   = params.Cg;
Cgsq = Cg^2;   
Hdot    = params.Hdot;
Roi  = params.Roi;
Beta   = params.Beta;
if (Roi~=0), 
    Ld = Cg/Roi;  % Nondimen. deformation scale {sqrt(g*H_0)/(f_0*L)}
end

shift = true;  % flag parameter for dif() and avg() to circshift result

% Coriolis parameter
if (Beta==0)
    FCOR.u = Roi;
    FCOR.v = Roi;
else
    yu = linspace(-L/2+DY/2,L/2-DY/2,ny-1);
    yv = linspace(-L/2,L/2-DY,ny-1);
    FCOR.u = meshgrid(Roi + Beta*yu,1:nx-1);      % f on u grid
    FCOR.v = meshgrid(Roi + Beta*yv,1:nx-1);      % f on v grid
end

% Do this for array F below, since this is how you shift to
% different grids in the RHS function
%
% Velocity fields
if (size(F_in.u,1)==nx&size(F_in.u,2)==ny-1)
    u = F_in.u;
elseif (params.geovel)  % u = -(1/f0)*dh/dy
    disp('Setting u=-(1/f)*dh/dy')
    u = zeros(nx,ny-1);
    dhy = zeros(nx-1,ny-1);
    dhy(:,2:ny-2) = (h(:,3:ny-1)-h(:,1:ny-3))/(2*DY);  % Get dh/dy on h points (i,j)
    dhy(:,1) = (h(:,2)-h(:,1))/DY;
    dhy(:,ny-1) = (h(:,ny-1)-h(:,ny-2))/DY;
    u = (1/Roi)*avg(dhy,1,PERIODIC.x,1);  % averaged onto u points
else
    u = zeros(nx,ny-1);
    disp('Setting u=0')
end
if (size(F_in.v,1)==nx-1&size(F_in.v,2)==ny)
    v = F_in.v;
elseif (params.geovel)
    disp('Setting v=(1/f0)*dh/dx')
    v = zeros(nx-1,ny);
    dhx = zeros(nx-1,ny-1);
    dhx(2:nx-2,:) = (h(3:nx-1,:)-h(1:nx-3,:))/(2*DX);  % Get dh/dx on h points 
    dhx(1,:) = (h(2,:)-h(1,:))/DX;
    dhx(nx-1,:) = (h(nx-1,:)-h(nx-2,:))/DX;
    v = (1/Roi)*avg(dhx,2,PERIODIC.y,1);
else
    v = zeros(nx-1,ny);
    disp('Setting v=0')
end

% Adaptive timestep and viscosity tuning factors
% umax should use actual wave speed, not Cg
umax = max([max(max(abs(u))) max(max(abs(v))) Cg]);
dt   = params.dttune*dr/umax;
nu   = params.Nu*dr^2/dt;
drag = params.Drag;

% Constant coefficients for RK3 timestepping
c1 = 1/3; c2 = 5/9; c3 = 15/16; c4 = 153/128; c5 = 8/15;

% Print some messages
disp('-------------')
params
if (Roi~=0), msg('Ld         :',Ld), end
msg('frame      :',frame)
msg('time       :',time)
msg('initial nu :',nu)
msg('initial dt :',dt)
disp('-------------')

% Store fields in single array for use in RK3
F = zeros(nx,ny,3);
R = F;       % RHS array
R1 = F;
F1 = F;
F2 = F;
zeta = [];
q = [];

F(:,1:end-1,1)       = u;
F(1:end-1,:,2)       = v;
F(1:end-1,1:end-1,3) = h-hb;

ke = h.*(sum(sum(u.^2))+sum(sum(v.^2)))/2; %% need an h multiplying that
ape = Cgsq*sum(sum((h).^2))/2; 
htot = sum(sum(h-hb));
t = time;

for n = 1:params.nt
  
    % Get rhs terms here so that zeta and q can be saved.
    R = rhs(F);
    
    % Save output fields and print diagnostics
    if ((mod(n,params.diagstep)==0)||(n==1)) 
        u = F(:,1:end-1,1);
        v = F(1:end-1,:,2);
        h = F(1:end-1,1:end-1,3)+hb;  
        time = [time t];
        ken = (sum(sum(u.^2))+sum(sum(v.^2)))/2;  ke = [ke ken];  % Append
        apen = Cgsq*sum(sum(h.^2))/2;             ape = [ape apen];
        htotn = sum(sum(h-hb));                   htot = [htot htotn];
        if (isnan(ken)|isinf(ken)|ken>.5/eps|~isreal(ken))
            msg('Blow up',ken)
            break
        end
        msg('')
        msg('Time     = ',t)
        msg('dt       = ',dt)
        msg('nu       = ',nu)
        msg('Total E  = ',ken+apen)
        %msg('Total H  = ',htotn)
    end
    if ((mod(n,params.savestep)==0)||(n==1))
        frame = frame+1;
        u = F(:,1:end-1,1);
        v = F(1:end-1,:,2);
        h = F(1:end-1,1:end-1,3)+hb;  
        q = (FCOR.v + zeta)./avg(avg(H,1,params.periodx,shift),2,params.periody,shift);

        if (params.writetofiles)
            write_field(h,hname,frame);
            write_field(q,qname,frame);
            write_field(zeta,zetaname,frame);
            write_field(u,uname,frame);
            write_field(v,vname,frame);
        else
            F_out(frame).u = u;
            F_out(frame).v = v;
            F_out(frame).h = h;    
            F_out(frame).q = q;    
            F_out(frame).zeta = zeta;    
        end
        msg('Wrote frame :',frame)
    end
    
    % RK3 timestepping
    R  = dt*R;
    F1 = F  + c1*R;
    R1 = dt*rhs(F1) - c2*R;
    F2 = F1 + c3*R1;
    F  = F2 + c5*(dt*rhs(F2) - c4*R1);
    
    % Make sure no normal flow BC enforced if not periodic
    if (~params.periodx)
        F(1,:,1)    = 0; % u = 0 at west
        F(end,:,1)  = 0; % u = 0 at east
    end
    if (~params.periody)
        F(:,1,2)    = 0; % v = 0 at south
        F(:,end,2)  = 0; % v = 0 at north
    end
    
    % Update clock
    t = t + dt;
    
    % Adapt dt and nu
    umax = max([max(max(abs(u))) max(max(abs(v))) Cg]);
    dt = params.dttune*dr/umax;
    nu = params.Nu*dr^2/dt;
    
end

% Save output
Diag_out.time = time;
Diag_out.ke = ke;
Diag_out.ape = ape;
Diag_out.htot = htot;

Params_out = params;
Params_out.nu = nu;
Params_out.dt = dt;

F_out(end).frame = frame;
F_out(end).time = time(end);
F_out(end).nx = nx;
F_out(end).ny = ny;
F_out(end).dx = DX;
F_out(end).dy = DY;
if (params.writetofiles)  %... then save restart to output structure F_out
    F_out.u = F(:,1:end-1,1);
    F_out.v = F(1:end-1,:,2);
    F_out.h = F(1:end-1,1:end-1,3);
    F_out.q = q;
    F_out.zeta = zeta;
end 


%-------------------------------------------------------------------
% Internal functions:  The 'end' at the end of the file means that
% all internal functions see variables from main, but not visa-versa
%-------------------------------------------------------------------

    function [Rlocal] = rhs(Flocal)

    persistent Rlocal Flocal B H % allocate memory and hold it
    
    % u v h and H are all nx-1 by ny-1 here.  Okay for periodic BC.
    % For walls, assume free-slip, so vorticity on walls is zero anyway.
    u = Flocal(1:end-1,1:end-1,1);
    v = Flocal(1:end-1,1:end-1,2);
    H = Flocal(1:end-1,1:end-1,3);      % h-hb;  Advect this in continuity eqn.
    h = H + hb;                    % Use this to calculate gradients for u,v
    
    % Vorticity on vorticity points 
    zeta = dif(v,1,params.periodx,shift)/DX - dif(u,2,params.periody,shift)/DY;
    
    % Make free slip BC in nonperiodic case
    if (~params.periodx), zeta(1,:) = 0; u(1,:) = 0; end
    if (~params.periody), zeta(:,1) = 0; v(:,1) = 0; end
    
    % Bernouli on h points
    B = Cgsq*h + ((avg(u,1,params.periodx)).^2 + (avg(v,2,params.periody)).^2)/2;
    
    % u_t = (f+zeta)*v - B_x + nu del2 u - drag*u
    
    %v_avxy = avg(avg(v,1,params.periodx,shift),2,params.periody);
    %B_x = dif(B,1,params.periodx,shift)/DX;
    %R(1:end-1,1:end-1,1) = v_avxy.*(FCOR.u + avg(zeta,2,params.periody)) - B_x + nulapu - r*u;
    %if (nu~=0), nutemp = nu*laplacian(u,DX,DY,params.periodx,params.periody); end
    
    Rlocal(1:end-1,1:end-1,1) = ...
        avg(avg(v,1,params.periodx,shift),2,params.periody).*(FCOR.u+avg(zeta,2,params.periody)) ...
        - dif(B,1,params.periodx,shift)/DX ...
        + nu*laplacian(u,DX,DY,params.periodx,params.periody) ...
        - drag*u;
    
    % v_t = -(f+zeta)*u - B_y + nu del2 v - drag*v
    
    %u_avxy = avg(avg(u,1,params.periodx),2,params.periody,shift);
    %B_y = dif(B,2,params.periody,shift)/DY;
    %R(1:end-1,1:end-1,2) = - u_avxy.*(FCOR.v + avg(zeta,1,params.periodx)) - B_y + nulapv - r*v;
    %if (nu~=0),  nutemp = nu*laplacian(v,DX,DY,params.periodx,params.periody); end
    
    Rlocal(1:end-1,1:end-1,2) = ...
        - avg(avg(u,1,params.periodx),2,params.periody,shift).*(FCOR.v+avg(zeta,1,params.periodx)) ...
        - dif(B,2,params.periody,shift)/DY ...
        + nu*laplacian(v,DX,DY,params.periodx,params.periody) ...
        - drag*v;
    
    
    %  h_t = -d(u(h-hb))/dx - d(v(h-hb))/dy + Hdot
    
    %H_avx = avg(H,1,params.periodx,shift);  % H averaged onto u points
    %H_avy = avg(H,2,params.periody,shift);  % H averaged onto v points
    %duHx = dif(u.*H_avx,1,params.periodx)/DX;
    %dvHy = dif(v.*H_avy,2,params.periody)/DY;
    %R(1:end-1,1:end-1,3) = -duHx - dvHy + Hdot;
    
    Rlocal(1:end-1,1:end-1,3) = ...
        -dif(u.*avg(H,1,params.periodx,shift),1,params.periodx)/DX ...
        -dif(v.*avg(H,2,params.periody,shift),2,params.periody)/DY + Hdot;
    
end

%---------------------------------------------------------------------
end % End of main function

