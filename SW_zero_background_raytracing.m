image_path = "images/";
write_file = false;

f = 3
gH = 1
Cg = sqrt(gH);

L = 2*pi;
nx = 512;

X = linspace(-L/2, L/2, nx);
[XX, YY] = meshgrid(X);

%streamfunction = @(x, y, t) (0.05*(x.^2 + y.^2));
%scheme = DifferenceScheme(streamfunction);

load ./ray_trace_sw/wavevort_231058_restart_frame100
scheme = SpectralScheme(Cg, f, L, nx, S);

Nparticles = 30;

x = zeros(Nparticles, 2);
k = zeros(Nparticles, 2);
t = 0;
for i=1:Nparticles
   k(i, :) = 3 * [cos(2*pi*i/Nparticles), sin(2*pi*i/Nparticles)]; 
   x(i, :) = L*(rand(1, 2)) - L/2;
end

U = scheme.U([XX, YY]);
speed = vecnorm(U, 2, 2);
U0 = max(speed(:));
C0 = sqrt(gH);
gH = Cg^2;
Fr = U0/C0
dx = L/nx;
dt = 0.5*dx/max(C0, U0);

Tend = 0.5/(f*Fr^2);
Nsteps = floor(Tend/dt)
% figure
% subplot(1,1,1);
% omega_s = spectrum(x, k, t);
% histogram(omega_s);
% figure
Omega_0 = omega(k, f, gH) + dot(scheme.U(x, t), k, 2);
error = zeros(Nsteps + 1, Nparticles);
solver_error = zeros(Nsteps + 1, Nparticles);

t_hist = zeros(Nsteps + 1, 1);
x_hist = zeros(Nsteps + 1, Nparticles, 2);
k_hist = zeros(Nsteps + 1, Nparticles, 2);
solver_x = zeros(Nsteps + 1, Nparticles, 2);
solver_k = zeros(Nsteps + 1, Nparticles, 2);

t_hist(1) = 0;
x_hist(1,:,:) = x;
k_hist(1,:,:) = k;

y0 = [x(:,1); x(:,2); k(:,1); k(:,2)];

rayfun = initialize_raytracing(scheme, f, gH, Nparticles);
tic
t_hist = dt * (0:Nsteps);
[t_hist, solver_y] = ode15s(rayfun, t_hist, y0);
toc
solver_x(:,:,1) = solver_y(:,1:Nparticles);
solver_x(:,:,2) = solver_y(:,Nparticles+1:2*Nparticles);
solver_k(:,:,1) = solver_y(:,2*Nparticles+1:3*Nparticles);
solver_k(:,:,2) = solver_y(:,3*Nparticles+1:4*Nparticles);

for i = 1:Nsteps
    Omega_abs = omega(squeeze(solver_k(i,:,:)), f, gH) + dot(scheme.U(squeeze(solver_x(i,:,:)), t), squeeze(solver_k(i,:,:)), 2);
    solver_error(i+1,:) = (Omega_abs - Omega_0) ./ Omega_0;
end

tic
fprintf("Simulation progress:  0.00%%")
for i = 1:Nsteps
   x_new = x + rk4(dt, k, x, t, @(k, x, t) (scheme.U(x, t) + grad_omega(k, f, gH)));
   k = k + rk4(dt, k, x, t, @(k, x, t) -scheme.grad_U_times_k(x, k, t));
   x = x_new;
   t = t + dt;
   t_hist(i+1) = t;
   x_hist(i+1,:,:) = x;
   k_hist(i+1,:,:) = k;
   Omega_abs = omega(k, f, gH) + dot(scheme.U(x, t), k, 2);
   error(i+1,:) = (Omega_abs - Omega_0) ./ Omega_0;
   if mod(i, 251) == 0
        fprintf("\b\b\b\b\b\b\b% 6.2f%%", i/Nsteps*100)
   end
end
toc


fmt = "yyyy-MM-dd HH_mm_ss";
filename = image_path + "QG_" + string(datetime("now"), fmt) + ".gif";

figure(1);
subplot(2, 1, 1);
p = scatter(x(:,1), x(:,2), 'k.');
axis image
hold on
p_solver = scatter(solver_x(1,:,1), solver_x(1,:,2), 'r.');
q = quiver(x(:,1), x(:,2), k(:,1), k(:,2), 0.5, 'k');
contour(XX, YY, scheme.streamfunction(XX, YY));
xlim([-L/2, L/2]);
ylim([-L/2, L/2]);
%legend("wavepacket")
xlabel("x");
ylabel("y");
grid on

subplot(2, 2, 3);
KK = scatter(k(:,1), k(:,2), 'k.');
xlabel("k");
ylabel("l");
grid on
axis image
xlim([-5, 5]);
ylim([-5, 5]);

subplot(2,2,4);
err_plot = plot(t_hist(1) * (f*Fr^2), error(1));
xlim([0, Tend * (f*Fr^2)]);
ylim([1.1*min(min(error)), 1.1*max(max(error))]);
%freq = sort(omega(squeeze(k_hist(1,:,:)), f, gH));
%W = plot(freq, energy(freq));

if(write_file)
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',inf);
end

anim = 0;
if(anim)
for i=1:40:Nsteps
       set(p, {'XData','YData'},{mod(x_hist(i,:,1) + L/2, L) - L/2, mod(x_hist(i,:,2) + L/2, L) - L/2});
       set(p_solver, {'XData', 'YData'},{mod(solver_x(i,:,1) + L/2, L) - L/2, mod(solver_x(i,:,2) + L/2, L) - L/2});
       set(q, {'XData', 'YData', 'UData', 'VData'}, {mod(x_hist(i,:,1) + L/2, L) - L/2, mod(x_hist(i,:,2) + L/2, L) - L/2, k_hist(i,:,1), k_hist(i,:,2)});
       set(KK, {'XData','YData'}, {k_hist(i,:,1), k_hist(i,:,2)});
       %freq = sort(omega(squeeze(k_hist(i,:,:)), f, gH));
       %set(W, {'XData', 'YData'}, {freq, energy(freq)});
       set(err_plot, {'XData','YData'}, {t_hist(1:i) * (f*Fr^2), error(1:i)})
       test = diff(error);
       
       %subplot(2, 2, 2);
       %subplot(1, 2, 2);
       %histogram(omega(k), 'Normalization', 'probability');
       %histogram(w.*sum(w' > w-5 & w' < w+5)');
       %scatter(omega(k)/f, energy(omega(k)));
       %xlim([0.1, 100]);
       %xlabel('\omega/f');
       %ylabel('e(\omega)');
       %set(gca,'xscale','log')
       %set(gca,'yscale','log')
       if(write_file)
           frame = getframe(1);
           im = frame2im(frame);
           [imind,cm] = rgb2ind(im,256);
           imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append');
       end
       pause(1/60);
end
end

% Plot errors in absolute magnitude
figure()
%rk4_error_plot = plot(t_hist * (f*Fr^2), error, "k");
title("Error in absolute frequency");
xlabel("t (1/(f*Fr^2)");
ylabel("\Delta\omega_a/\omega_0");
hold on
solver_error_plot = plot(t_hist * (f*Fr^2), solver_error);
%legend([rk4_error_plot(1), solver_error_plot(1)], "RK4 error", "MATLAB solver error");
xlabel("t (1/(f*Fr^2)");
ylabel("\Delta\omega_a/\omega_0");

figure()
title("Accumulation plot")
plot(k_hist(:,:,1), k_hist(:,:,2))
hold on
plot(k_hist(1,:,1), k_hist(1,:,2), 'k', 'LineWidth', 2.5)
hold off
xlabel("k_1");
ylabel("k_2");

function rayode = initialize_raytracing(scheme, f, gH, Nparticles)
    function dydt = odefun(t,y)
       x = [y(1:Nparticles), y(Nparticles+1:2*Nparticles)];
       k = [y(2*Nparticles+1:3*Nparticles), y(3*Nparticles+1:4*Nparticles)];
       dxdt = scheme.U(x, t) + grad_omega(k, f, gH);
       dkdt = -scheme.grad_U_times_k(x, k, t);
       dydt = zeros(4*Nparticles, 1);
       dydt(0*Nparticles + 1:1*Nparticles) = dxdt(:, 1);
       dydt(1*Nparticles + 1:2*Nparticles) = dxdt(:, 2);
       dydt(2*Nparticles + 1:3*Nparticles) = dkdt(:, 1);
       dydt(3*Nparticles + 1:4*Nparticles) = dkdt(:, 2);
       
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