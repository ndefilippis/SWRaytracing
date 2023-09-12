addpath ./qg_flow_ray_trace

image_path = "images/";
write_file = false;

run_log_file = "analysis/job-36976465/run-4/run.log";
[nx, ~, f, Cg, Ug] = parse_data(run_log_file);
K_d2 = f/Cg;

kmax = nx/2-1;
[kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
K2 = kx_.^2 + ky_.^2;

rng(123)
L = 2*pi;
dx = L/nx;

X = linspace(-L/2, L/2, nx);
num_interp_points = 4;
[XX, YY] = meshgrid(X);

%streamfunction = @(x, y, t) (0.5*(x.^2 + y.^2 - x.^4));
%scheme = DifferenceScheme(streamfunction);

times = read_field("analysis/pv_time");
q = read_field("analysis/pv", nx, nx, 1, [2000]);
background_flow = grid_U(g2k(q), K_d2, K2, kx_, ky_);
%S = zeros(nx, nx, 3);
%S(:,:,3) = 0.05;
scheme = SpectralScheme(L, nx, k2g(-g2k(q)./(K_d2 + K2)));

Nparticles = 10;

x = zeros(1, 2, Nparticles);
k = zeros(1, 2, Nparticles);
t = 0;
for i=1:Nparticles
   k(1, :, i) = 3 * [cos(2*pi*i/Nparticles), sin(2*pi*i/Nparticles)]; 
   x(1, :, i) = L*(rand(1, 2)) - L/2;
end

U = scheme.U([XX(:), YY(:)]);
speed = vecnorm(U, 2, 2);
dx = L/nx;
U0 = max(speed(:));
gH = Cg^2;
Fr = U0/Cg
dt = 0.1*dx/max(Cg, U0);

Tend = 1/(f*Fr^2);
Nsteps = floor(Tend/dt)
% figure
% subplot(1,1,1);
% omega_s = spectrum(x, k, t);
% histogram(omega_s);
% figure
Omega_0 = squeeze(omega(k, f, gH) + dot(scheme.U(x, t), k, 2));

t_hist = zeros(Nsteps + 1, 1);
x_hist = zeros(Nsteps + 1, 2, Nparticles);
k_hist = zeros(Nsteps + 1, 2, Nparticles);
solver_x = zeros(Nsteps + 1, 2, Nparticles);
solver_k = zeros(Nsteps + 1, 2, Nparticles);

t_hist(1) = 0;
x_hist(1,:,:) = x;
k_hist(1,:,:) = k;

y0 = [squeeze(x)' squeeze(k)'];

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-7);
method = @ode23;

rayfun = initialize_raytracing(scheme, f, gH, Nparticles);
t_hist = dt * (0:Nsteps);
fprintf(func2str(method)+" with 1e-5 tol:\n")
tic
[t_hist, solver_y] = method(rayfun, t_hist, y0(:), opts);
toc
solver_x(:,1,:) = solver_y(:,1:Nparticles);
solver_x(:,2,:) = solver_y(:,Nparticles+1:2*Nparticles);
solver_k(:,1,:) = solver_y(:,2*Nparticles+1:3*Nparticles);
solver_k(:,2,:) = solver_y(:,3*Nparticles+1:4*Nparticles);

w = squeeze(omega(solver_k, f, gH));
Omega_abs = squeeze(omega(solver_k, f, gH) + dot(scheme.U(solver_x, t), solver_k, 2));
solver_error = (Omega_abs' - Omega_0) ./ Omega_0;

gap = 500;
samples = w(1000:gap:end,:);
samples = samples(:);
hist(samples);


% for i=1:4:Nsteps
%   subplot(2, 1, 1);
%   contour(XX, YY, scheme.streamfunction(XX, YY));
%   axis image
%   hold on
%   scatter(mod(solver_x(i,1,:) + L/2, L) - L/2, mod(solver_x(i,2,:) + L/2, L) - L/2, 30, 'k.');
%   hold off
%   subplot(2, 1, 2);
%   scatter(solver_k(i,1,:), solver_k(i,2,:));
%   axis image;
%   pause(1/10); 
% end


%tic
%fprintf("Simulation progress:  0.00%%")
%for i = 1:Nsteps
%   x_new = x + rk4(dt, k, x, t, @(k, x, t) (scheme.U(x, t) + grad_omega(k, f, gH)));
%   k = k + rk4(dt, k, x, t, @(k, x, t) -scheme.grad_U_times_k(x, k, t));
%   x = x_new;
%   t = t + dt;
%   t_hist(i+1) = t;
%   x_hist(i+1,:,:) = x;
%   k_hist(i+1,:,:) = k;
%   Omega_abs = omega(k, f, gH) + dot(scheme.U(x, t), k, 2);
%   error(i+1,:) = (Omega_abs - Omega_0) ./ Omega_0;
%   if mod(i, 251) == 0
%        fprintf("\b\b\b\b\b\b\b% 6.2f%%", i/Nsteps*100)
%   end
%end
%toc

% Plot errors in absolute magnitude
figure()
solver_error_plot = plot(t_hist * (f*Fr^2), solver_error, "k");
title("Error in absolute frequency");
xlabel("t (1/(f*Fr^2)");
ylabel("\Delta\omega_a/\omega_0");

function rayode = initialize_raytracing(scheme, f, gH, Nparticles)
    function dydt = odefun(t,y)
       y = reshape(y,Nparticles,[]);
       x = y(:,1:2);
       k = y(:,3:4);
       dxdt = scheme.U(x, t) + grad_omega(k, f, gH);
       dkdt = -scheme.grad_U_times_k(x, k, t);
       dydt = [dxdt dkdt];
       dydt = dydt(:);       
    end
    rayode = @odefun;
end

function [resolution, Npackets, f, Cg, Ug] = parse_data(filename) 
    [fid,status]=fopen(filename);
    if(status)
        disp(status);
        return;
    end
    resolution_cell = textscan(fid,'Resolution: %fx%f',1,'delimiter','\n', 'headerlines', 10);    
    Npackets_cell = textscan(fid, "Number of packets: %d",1);
    f_cell = textscan(fid, "Coriolis parameter: %f",'delimiter','\n', 'headerlines', 7);
    Cg_cell = textscan(fid,  "Group velocity: %f", 1);
    Ug_cell = textscan(fid, "Background velocity (parameter,computed): (%f,%f)", 1);
    resolution = resolution_cell{1};
    Npackets = Npackets_cell{1};
    f = f_cell{1};
    Cg = Cg_cell{1};
    Ug = Ug_cell{1};
end

function omega_s=spectrum(x, k, t)
    [psi_xx, psi_xy, psi_yy] = grad_U(x, k, t);
    D = psi_xy.^2 - psi_xx.*psi_yy;
    D = D(:);
    
    eig_k = ones(size(D));
    eig_l = -(-psi_xy - sqrt(D)).*eig_k ./ psi_xx;
    eig = [eig_k, eig_l];
    eig = 30 .* eig ./ (1 + eig_l.^2);
    eig = eig(D > 0);
    omega_s = omega(eig);
end

function w=omega(k, f, gH)
    w = sqrt(f*f + gH.*dot(k, k, 2));
end

function cg=grad_omega(k, f, gH)
    cg = gH * k ./ sqrt(f*f + gH.*(dot(k, k, 2)));
end

function update=rk4(dt, k, x, t, f)
    k1 = f(k, x, t);
    k2 = f(k + dt * k1/2, x + dt * k1/2, t + dt/2);
    k3 = f(k + dt * k2/2, x + dt * k2/2, t + dt/2);
    k4 = f(k + dt * k3, x + dt * k3, t + dt);
    
    update = dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

function e=energy(omega)
    e = omega.*sum(omega' > omega-3 & omega' < omega+3)';
end
