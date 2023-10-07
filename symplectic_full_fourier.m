addpath ./qg_flow_ray_trace

run_log_file = "analysis/job-36976465/run-4/run.log";
[nx, ~, f, Cg, Ug] = parse_data(run_log_file);
K_d2 = f/Cg;

kmax = nx/2-1;
[kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
K2 = kx_.^2 + ky_.^2;

rng(123)
L = 2*pi;
dx = L/nx;

X = linspace(0, L, nx);
[XX, YY] = meshgrid(X);

q = read_field("analysis/pv", nx, nx, 1, [2000]);

scheme = SpectralScheme(L, nx, k2g(-g2k(q)./(K_d2 + K2)));

Nparticles = 10;
x = zeros(1, 2, Nparticles);
k = zeros(1, 2, Nparticles);
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

Tend = 10/(f*Fr^2);


Omega_0 = omega(k, f, gH) + dot(scheme.U(x, 0), k, 2);

tic
[solver_x, solver_k, solver_t] = ode_symplectic(x, k, dt, Tend, f, gH, scheme);
toc

% for i=1:10:length(solver_t)
%    k = squeeze(solver_k(i, 1, :));
%    l = squeeze(solver_k(i, 2, :));
%    scatter(k, l, 20, 'k.');
%    pause(1/30);
% end

w = squeeze(omega(solver_k, f, gH));
Omega_abs = omega(solver_k, f, gH) + dot(scheme.U(solver_x, 0), solver_k, 2);
solver_error = squeeze((Omega_abs - Omega_0) ./ Omega_0);
plot(solver_t * (f*Fr^2), solver_error, 'k');
title("Error in absolute frequency");
xlabel("t (1/(f*Fr^2)");
ylabel("\Delta\omega_a/\omega_0");

function w=omega(k, f, gH)
    w = sqrt(f*f + gH.*dot(k, k, 2));
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