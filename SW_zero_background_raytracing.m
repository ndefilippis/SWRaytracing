image_path = "images/";
write_file = false;

rng(123)

f = 3
gH = 1
Cg = sqrt(gH);

L = 2*pi;
nx = 512;

X = linspace(-L/2, L/2, nx);
[XX, YY] = meshgrid(X);

%streamfunction = @(x, y, t) (0.05*(x.^2 + y.^2 - x.^4));
%scheme = DifferenceScheme(streamfunction);

load ./ray_trace_sw/wavevort_231058_restart_frame100
%S = zeros(nx, nx, 3);
%S(:,:,3) = 0.05;
scheme = SpectralScheme(Cg, f, L, nx, S);

Nparticles = 10;

x = zeros(Nparticles, 2);
k = zeros(Nparticles, 2);
t = 0;
for i=1:Nparticles
   k(i, :) = 3 * [cos(2*pi*i/Nparticles), sin(2*pi*i/Nparticles)]; 
   x(i, :) = L*(rand(1, 2)) - L/2;
end

U = scheme.U([XX(:), YY(:)]);
speed = vecnorm(U, 2, 2);
U0 = max(speed(:));
C0 = sqrt(gH);
gH = Cg^2;
Fr = U0/C0
dx = L/nx;
dt = 0.5*dx/max(C0, U0);

Tend = 1/(f*Fr^2);
Nsteps = floor(Tend/dt)
% figure
% subplot(1,1,1);
% omega_s = spectrum(x, k, t);
% histogram(omega_s);
% figure
Omega_0 = omega(k, f, gH) + dot(scheme.U(x, t), k, 2);
error = zeros(Nsteps + 1, Nparticles);
solver_error = zeros(Nsteps + 1, Nparticles);
solver2_error = zeros(Nsteps + 1, Nparticles);

t_hist = zeros(Nsteps + 1, 1);
x_hist = zeros(Nsteps + 1, Nparticles, 2);
k_hist = zeros(Nsteps + 1, Nparticles, 2);
solver_x = zeros(Nsteps + 1, Nparticles, 2);
solver_k = zeros(Nsteps + 1, Nparticles, 2);
solver_x2 = zeros(Nsteps + 1, Nparticles, 2);
solver_k2 = zeros(Nsteps + 1, Nparticles, 2);

t_hist(1) = 0;
x_hist(1,:,:) = x;
k_hist(1,:,:) = k;

y0 = [x(:,1); x(:,2); k(:,1); k(:,2)];


opts1 = odeset('RelTol', 1e-3, 'AbsTol', 1e-5);
opts2 = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
method1 = @ode23;
method2 = @ode23;

rayfun = initialize_raytracing(scheme, f, gH, Nparticles);
t_hist = dt * (0:Nsteps);
fprintf(func2str(method1)+" with 1e-3 tol:\n")
tic
[t_hist, solver_y] = method1(rayfun, t_hist, y0, opts1);
toc
fprintf(func2str(method2)+" with 1e-5 tol:\n")
tic
[t_hist, solver_y2] = method2(rayfun, t_hist, y0, opts);
toc
solver_x(:,:,1) = solver_y(:,1:Nparticles);
solver_x(:,:,2) = solver_y(:,Nparticles+1:2*Nparticles);
solver_k(:,:,1) = solver_y(:,2*Nparticles+1:3*Nparticles);
solver_k(:,:,2) = solver_y(:,3*Nparticles+1:4*Nparticles);
solver_x2(:,:,1) = solver_y2(:,1:Nparticles);
solver_x2(:,:,2) = solver_y2(:,Nparticles+1:2*Nparticles);
solver_k2(:,:,1) = solver_y2(:,2*Nparticles+1:3*Nparticles);
solver_k2(:,:,2) = solver_y2(:,3*Nparticles+1:4*Nparticles);

for i = 1:Nsteps
    Omega_abs = omega(squeeze(solver_k(i,:,:)), f, gH) + dot(scheme.U(squeeze(solver_x(i,:,:)), t), squeeze(solver_k(i,:,:)), 2);
    solver_error(i+1,:) = (Omega_abs - Omega_0) ./ Omega_0;
    Omega_abs = omega(squeeze(solver_k2(i,:,:)), f, gH) + dot(scheme.U(squeeze(solver_x2(i,:,:)), t), squeeze(solver_k2(i,:,:)), 2);
    solver2_error(i+1,:) = (Omega_abs - Omega_0) ./ Omega_0;
end

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
hold on
solver2_error_plot = plot(t_hist * (f*Fr^2), solver2_error, 'r--');
legend([solver_error_plot(1), solver2_error_plot(1)], func2str(method1) + " error", func2str(method2) + " error");
xlabel("t (1/(f*Fr^2)");
ylabel("\Delta\omega_a/\omega_0");

function S = sparsity(Npackets)
    S = zeros(4*Npackets);
    for i=1:Npackets
       for field=0:3
           for field2=0:3
            S(i + Npackets*field, i + Npackets*field2) = 1;
            S(i + Npackets*field2, i + field*Npackets) = 1;
           end
       end
    end
end

function rayode = initialize_raytracing(scheme, f, gH, Nparticles)
    function dydt = odefun(t,y)
       x = [y(1:Nparticles), y(Nparticles+1:2*Nparticles)];
       k = [y(2*Nparticles+1:3*Nparticles), y(3*Nparticles+1:4*Nparticles)];
       dxdt = scheme.U(x, t) + grad_omega(k, f, gH);
       dkdt = -scheme.grad_U_times_k(x, k, t);
       dydt = zeros(4*Nparticles, 1);
       dydt(0*Nparticles + 1:1*Nparticles,:) = dxdt(:, 1);
       dydt(1*Nparticles + 1:2*Nparticles,:) = dxdt(:, 2);
       dydt(2*Nparticles + 1:3*Nparticles,:) = dkdt(:, 1);
       dydt(3*Nparticles + 1:4*Nparticles,:) = dkdt(:, 2);
       
    end
    rayode = @odefun;
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
