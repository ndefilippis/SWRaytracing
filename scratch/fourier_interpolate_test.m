addpath ../rsw

L = 2*pi;
nx = 32;
x = linspace(0, L, nx);
[X, Y] = meshgrid(x, x);

n = 1;
N = 2*n + 1;
amp = rand(N, N)/(N^2);
phase = L*rand(N, N);

Nparticles = 10;
x0 = zeros(1, 2, Nparticles);
k0 = zeros(1, 2, Nparticles);
for i=1:Nparticles
   k0(1, :, i) = 3 * [cos(2*pi*i/Nparticles), sin(2*pi*i/Nparticles)]; 
   x0(1, :, i) = L*(rand(1, 2)) - L/2;
end

f = 3;
Cg = 1;
gH = Cg^2;
U = Velocity(X(:), Y(:), amp, phase, n);
speed = vecnorm(U, 2, 2);
dx = L/nx;
U0 = max(speed(:));
dt = 0.1*dx/max(Cg, U0);

Fr = U0/sqrt(gH)

Tend = 10/(f*Fr^2);

Omega_0 = omega(k0, f, gH) + dot(Velocity(x0(1, 1, :), x0(1, 2, :), amp, phase, n), k0(1, :, :), 2);

Nsteps = floor(Tend / dt);
solver_x = zeros([Nsteps, size(x0, [2,3])]);
solver_k = zeros([Nsteps, size(k0, [2,3])]);
solver_t = zeros(Nsteps, 1);
solver_x(1, :, :) = x0;
solver_k(1, :, :) = k0;
solver_t(1) = 0;

for i=2:Nsteps
    [x0, k0] = leapfrog_method(x0, k0, dt, @phi1, @phi2, f, gH, amp, phase, n);
    solver_x(i, :, :) = x0;
    solver_k(i, :, :) = k0;
    solver_t(i) = (i-1)*dt;
end

w = squeeze(omega(solver_k, f, gH));
Omega_abs = omega(solver_k, f, gH) + dot(Velocity(solver_x(:,1,:), solver_x(:,2,:), amp, phase, n), solver_k, 2);
solver_error = squeeze((Omega_abs - Omega_0) ./ Omega_0);
plot(solver_t * (f*Fr^2), solver_error, 'k');
title("Error in absolute frequency");
xlabel("t (1/(f*Fr^2)");
ylabel("\Delta\omega_a/\omega_0");

function [x, k] = leapfrog_method(x0, k0, dt, phi1, phi2, f, gH, amp, phase, n)
        [x1, k1] = phi1(x0, k0, dt/2, f, gH);
        [x2, k2] = phi2(x1, k1, dt, amp, phase, n);
        [x, k] = phi1(x2, k2, dt/2, f, gH);
end
    
function [x, k] = phi1(x0, k0, dt, f, gH)
    x = x0 + dt * group_velocity(k0, f, gH);
    k = k0;
end

function w = omega(k, f, gH)
    w = sqrt(f^2 + gH*dot(k, k, 2));
end

function cg = group_velocity(k, f, gH)
    cg = gH * k./omega(k, f, gH);
end

function [x, k] = phi2(x0, k0, dt, amp, phase, n)
    xx = x0(:,1,:);
    yy = x0(:,2,:);
    kk = k0(:,1,:);
    ll = k0(:,2,:);
    x = xx;
    y = yy;
    k = kk;
    l = ll;
    for K=-n:n
        for L=-n:n
            if K == 0 && L == 0
                continue
            end
            c = K*xx + L*yy;
            a = L*kk + K*ll;
            
            x = x - dt * L * real(1i * amp(K+n+1, L+n+1)*exp(1i*(c + phase(K+n+1,L+n+1))));
            y = y + dt * K * real(1i * amp(K+n+1, L+n+1)*exp(1i*(c + phase(K+n+1,L+n+1))));
            k = k - dt * K * real(-1 * amp(K+n+1, L+n+1)*exp(1i*(c + phase(K+n+1,L+n+1))) .* a);
            l = l - dt * L * real(-1 * amp(K+n+1, L+n+1)*exp(1i*(c + phase(K+n+1,L+n+1))) .* a);
        end
    end
    x = [x y];
    k = [k l];
end

function psi=streamfunction(arg, amp, phase, n)
    psi = 0 * arg;
    for k=-n:n
        for l=-n:n
            if k == 0 && l == 0
                continue
            end
            psi = psi + amp(k+n+1, l+n+1)*real(exp(1i*(arg + phase(k+n+1,l+n+1))));
        end
    end
end

function vel = Velocity(x, y, amp, phase, n)
    vel = [0 * x, 0 * y];
    for k=-n:n
        for l=-n:n
            if k == 0 && l == 0
                continue
            end
            vel = vel + amp(k+n+1, l+n+1)*real(exp(1i*(k*x + l*y + phase(k+n+1,l+n+1))) .* [-1i*l, 1i*k]);
        end
    end
end

