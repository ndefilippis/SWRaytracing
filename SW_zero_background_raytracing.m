load C:/Users/ndefi/Downloads/ray_trace_sw/wavevort_231058_restart_frame100
% Contains state variable S(:,:,:) with [u,v,eta]=S(:,:,1:3)
% Use method from wavevortdecomp.m to extract geostrophic mode
% velocities ug, vg, and compute gradients spectrally
addpath C:/Users/ndefi/Downloads/rsw/
addpath C:/Users/ndefi/Downloads/ray_trace_sw/

image_path = "images/";

L = 2*pi;
nx = 512;
X = linspace(-L/2, L/2, nx);
[XX, YY] = meshgrid(X);
% S(:,:,1) = 0;
% S(:,:,2) = 0;
% S(:,:,3) = (cos(XX).^2 + cos(YY).^2);

streamfunction = @(x, y, t) (x.^2 + y.^2);

%scheme = DifferenceScheme(streamfunction);
scheme = SpectralScheme(Cg, f, L, nx, S);

fmt = "yyyy-MM-dd HH_mm_ss";
filename = image_path + "QG_" + string(datetime("now"), fmt) + ".gif";
write_file = false;
Nparticles = 30;

x = zeros(Nparticles, 2);
k = zeros(Nparticles, 2);
t = 0;
for i=1:Nparticles
   k(i, :) = 30 * [cos(2*pi*i/Nparticles), sin(2*pi*i/Nparticles)]; 
   x(i, :) = L*(rand(1, 2)) - L/2;
end

U = scheme.U([XX, YY]);
speed = vecnorm(U, 2, 2);
U0 = max(speed(:));
C0 = Cg;
gH = Cg^2;
Fr = U0/C0;
dx = L/nx;
dt = 0.3*dx/max(C0, U0);

figure(1);
subplot(2, 1, 1);
p = scatter(x(:,1), x(:,2), 'k.');
axis image
hold on
q = quiver(x(:,1), x(:,2), k(:,1), k(:,2), 0.5, 'k');
contour(XX, YY, scheme.streamfunction(XX, YY));
xlim([-L/2, L/2]);
ylim([-L/2, L/2]);
legend("wavepacket")
xlabel("x");
ylabel("y");
grid on

subplot(2, 2, 3);
KK = scatter(k(:,1), k(:,2), 'k.');
xlabel("k");
ylabel("l");
grid on
axis image
xlim([-50, 50]);
ylim([-50, 50]);

if(write_file)
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',inf);
end

Tend = 1/(f*Fr^2);
Nsteps = floor(Tend/dt);
% figure
% subplot(1,1,1);
% omega_s = spectrum(x, k, t);
% histogram(omega_s);
% figure
Omega_0 = omega(k, f, gH) + dot(scheme.U(x, t), k, 2);
error = zeros(Nsteps + 1, Nparticles);
t_hist = zeros(Nsteps + 1, 1);
t_hist(1) = 0;
for i = 1:Nsteps
   x = x + rk4(dt, k, x, t, @(k, x, t) (scheme.U(x, t) + grad_omega(k, f, gH)));
   k = k + rk4(dt, k, x, t, @(k, x, t) -scheme.grad_U_times_k(x, k, t));
   t = t + dt;
   t_hist(i+1) = t;
   Omega_abs = omega(k, f, gH) + dot(scheme.U(x, t), k, 2);
   error(i+1,:) = (Omega_abs - Omega_0) ./ Omega_0;
   if mod(i, 40) == 1
       set(p, {'XData','YData'},{mod(x(:,1) + L/2, L) - L/2, mod(x(:,2) + L/2, L) - L/2});
       set(q, {'XData', 'YData', 'UData', 'VData'}, {mod(x(:,1) + L/2, L) - L/2, mod(x(:,2) + L/2, L) - L/2, k(:,1), k(:,2)});
       set(KK, {'XData','YData'}, {k(:,1), k(:,2)});
       
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
figure()
plot(t_hist * (f*Fr^2), error);
xlabel("t (1/(f*Fr^2)");
ylabel("\Delta\omega_a/\omega_0");


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
