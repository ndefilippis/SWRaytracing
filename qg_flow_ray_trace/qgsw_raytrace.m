function qgsw_raytrace(nx, Npackets, near_inertial_factor, T_days, packet_delay_days, U_g, f, Cg)
% Input:
% nx: resolution of QG flow
% Npackets: number of packets to advect
% near_inertial_factor: How close to f the initial wavenumbers are
% T_days: number of days to simulate
% packet_delay_days: Number of days to wait before simulation packets
% U_g: amplitude of background geostrophic flow
% f: Coriolis parameter
% Cg: Group velocity of the waves (equal to sqrt(gH))

% Set up domain
L = 2*pi;
dx = L/nx;
x = linspace(-L/2, L/2, nx);
[X, Y] = meshgrid(x, x);

kmax = nx/2-1;
[kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
K2 = kx_.^2 + ky_.^2;

% Simulation parameters
rng(146);
beta = 0;
K_d2 = f/Cg;
T = T_days/f;
CFL_fraction = 0.1;

% Output parameters
steps_per_save = 50; % Has to be bigger than 3 or it wont save the initial AB steps
packet_delay = packet_delay_days / f;
packet_steps_per_save = 5;
pv_filename = 'pv';
pv_time_filename = 'pv_time';
packet_x_filename = 'packet_x';
packet_k_filename = 'packet_k';
packet_time_filename = 'packet_time';

% Create log levels
LOG_ERROR = 0;
LOG_INFO = 1;
LOG_VERBOSE = 2;
log_message = create_logger(LOG_INFO);


% Set up initial conditions
t = 0;

q = initial_q(X, Y, U_g, K_d2);
qk = g2k(q);

packet_x = zeros(Npackets, 2);
packet_k = zeros(Npackets, 2);
for i=1:Npackets
   packet_k(i, :) = near_inertial_factor * f * [cos(2*pi*i/Npackets), sin(2*pi*i/Npackets)]; 
   packet_x(i, :) = L*(rand(1, 2)) - L/2;
end

% Compute time step and Froude number from inital PV flow
flow = grid_U(qk, K_d2, K2, kx_, ky_);
speed2 = flow.u.^2 + flow.v.^2;
U0 = sqrt(max(speed2(:)));
Fr = U0/Cg;

dt = CFL_fraction*dx/U0;

Nsteps = ceil(T/dt);
packet_step_start = ceil(packet_delay / dt);

% Write out parameters:
log_message("Resolution: %dx%d\n", LOG_INFO, nx, nx);
log_message("Number of packets: %d\n", LOG_INFO, Npackets);
log_message("Initial wavenumber radius: %f\n", LOG_INFO, near_inertial_factor * f);
log_message("Time step: %f\n", LOG_INFO, dt);
log_message("Simulation time: %f\n", LOG_INFO, T);
log_message("Spin-up time: %f\n", LOG_INFO, packet_delay);
log_message("Steps per save: %d\n", LOG_INFO, steps_per_save);
log_message("Steps per packet save: %d\n", LOG_INFO, packet_steps_per_save);
log_message("Coriolis parameter: %f\n", LOG_INFO, f);
log_message("Group velocity: %f\n", LOG_INFO, Cg);
log_message("Background velocity (parameter,computed): (%f,%f)\n", LOG_INFO, U_g, U0);
log_message("Froude Number: %f\n", LOG_INFO, Fr);
log_message("Deformation wavenumber: %f\n", LOG_INFO, K_d2);


t_background_save = zeros(1 + floor(Nsteps / steps_per_save), 1);
q_save = zeros(nx, nx, 1 + floor(Nsteps / steps_per_save));
x_save = zeros(Npackets, 2, 1 + floor((Nsteps - packet_step_start + 1) / packet_steps_per_save));
k_save = zeros(Npackets, 2, 1 + floor((Nsteps - packet_step_start + 1) / packet_steps_per_save));
t_packet_save = zeros(1 + floor((Nsteps - packet_step_start + 1) / packet_steps_per_save), 1);

frame = 1;
packet_frame = 1;
q_save(:, :, frame) = q;
x_save(:, :, packet_frame) = packet_x;
k_save(:, :, packet_frame) = packet_k;

% Write initial positions
write_field(packet_x, packet_x_filename, packet_frame);
write_field(packet_k, packet_k_filename, packet_frame);
write_field(dt * (packet_step_start - 1), packet_time_filename, packet_frame);

write_field(q, pv_filename, frame);
write_field(t, pv_time_filename, frame);

AB_order = 3;
Qn_minus = zeros(2*kmax+1, kmax+1, AB_order - 1);

tic
Ef = filter(kx_, ky_, dx);
% Initialize Adams-Bashforth inital steps using:
% 1. Forward Euler
% 2. Second order AB
log_message("Simulation progress:  0.00%%", LOG_VERBOSE)
for step=1:Nsteps
   prev_qk = qk;
   if(step == 1)
       Qn = update(qk, K2, K_d2, beta, kx_, ky_);
       dq = dt*Qn;
   elseif(step == 2)
       Qn = update(qk, K2, K_d2, beta, kx_, ky_);
       dq = dt/2*(3*Qn - Qn_minus(:,:,1));
   else
       Qn = update(qk, K2, K_d2, beta, kx_, ky_);
       dq = dt/12*(23*Qn - 16*Qn_minus(:,:,1) + 5*Qn_minus(:,:,2));
   end
   t = t + dt;
   Qn_minus(:,:,2) = Qn_minus(:,:,1);
   Qn_minus(:,:,1) = Qn;
   qk = qk + dq;
   qk = Ef .* qk;% At some point, try to get hyperdiffusion working here instead/also
   
   % Do wavepacket advection
   if(t > packet_delay)
       background_flow1 = grid_U(prev_qk, K_d2, K2, kx_, ky_);
       background_flow2 = grid_U(qk, K_d2, K2, kx_, ky_);
       ray_ode = generate_raytracing_ode(background_flow1, background_flow2, Npackets, f, Cg, dt, dx);
       y0 = ode_xk2y(Npackets, packet_x, packet_k);
       % opts = odeset('Jpattern', S);
       % Changed from ode15s to od23 since we got better results for the
       % static case. This means we also don't need to compute or use the
       % Jacobian
       [~, solver_y] = ode23(ray_ode, [0, dt], y0);
       [packet_x, packet_k] = ode_y2xk(Npackets, solver_y(end,:)');
   end
   
   if(t > packet_delay && mod(step - packet_step_start + 1, packet_steps_per_save) == 0)
      packet_frame = packet_frame + 1; 
      x_save(:, :, packet_frame) = mod(packet_x + L/2, L) - L/2;
      k_save(:, :, packet_frame) = packet_k;
      t_packet_save(packet_frame) = t;
      
      % Write fields
      write_field(mod(packet_x + L/2, L) - L/2, packet_x_filename, packet_frame);
      write_field(packet_k, packet_k_filename, packet_frame);
      write_field(t, packet_time_filename, packet_frame);
   end
   
   if(mod(step, steps_per_save) == 0)
       frame = frame + 1;
       q = k2g(qk);
       q_save(:,:,frame) = q;
       t_background_save(frame) = t;
       write_field(q, pv_filename, frame);
       write_field(t, pv_time_filename, frame);
   end
   if mod(step, 51) == 0
        log_message("\b\b\b\b\b\b\b% 6.2f%%", LOG_VERBOSE, step/Nsteps*100)
   end
end
log_message("\b\b\b\b\b\b\b100.00%%\n", LOG_VERBOSE);
time_elapsed = toc;
log_message("Real time elapsed: %.3f seconds\n", LOG_INFO, time_elapsed);
end

function log_message = create_logger(max_log_level)
    function log_func(output_string, log_level, varargin)
        if(log_level <= max_log_level)
            fprintf(output_string, varargin{:}); 
        end
    end
    log_message = @log_func;
end

function q=initial_q(X, Y, a_g, K_d2)
    % Inital background PV
    % set as a ring of intermediate wavenumbers
    k_min = 1;
    k_max = 3;
    q = 0*X;
    U = 0*X;
    V = 0*X;
    phase = 2*pi*rand(2*k_max + 1, 2*k_max + 1);
    for k=-k_max:k_max
        for l=-k_max:k_max
            if k_min^2 < k^2 + l^2 <= k_max^2
                % Keep track of total eta to normalize after the face
                wave_phase = k*X + l*Y + phase(k+k_max + 1, l+k_max + 1);
                % The constant here is to satisfy geostrophy
                U = U - l*sin(wave_phase);
                V = V + k*sin(wave_phase);
                q = q - (K_d2 + k^2 + l^2)*cos(wave_phase);
            end
        end
    end
    speed2 = U.^2 + V.^2;
    q = a_g/sqrt(max(speed2(:))) * q;
end

function Ef=filter(kx_, ky_, dx)
    Ef = ones(size(kx_));
    kstar = sqrt((kx_*dx).^2 + (ky_*dx).^2);
    
    kc = 0.75*pi;
    const = log(1e-15)/(0.25*pi)^4;
    result = exp(const * (kstar - kc).^4);
    Ef(kstar >= kc) = result(kstar >= kc);
end

% Helper functions to turn x-k packet variables to a single column variable
function [x, k] = ode_y2xk(Npackets, y)
    x = [y(0*Npackets+1:1*Npackets), y(1*Npackets+1:2*Npackets)];
    k = [y(2*Npackets+1:3*Npackets), y(3*Npackets+1:4*Npackets)];
end

function y = ode_xk2y(Npackets, x, k)
    y = zeros(4*Npackets, 1);
    y(0*Npackets + 1:1*Npackets) = x(:, 1);
    y(1*Npackets + 1:2*Npackets) = x(:, 2);
    y(2*Npackets + 1:3*Npackets) = k(:, 1);
    y(3*Npackets + 1:4*Npackets) = k(:, 2);
end

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

function rayode=generate_raytracing_ode(background_flow1, background_flow2, Npackets, f, Cg, tmax, h)
    function dydt=odefun(t, y)
       [x, k] = ode_y2xk(Npackets, y);
       [U, nablaU] = interpolate_U(background_flow1, background_flow2, t/tmax, x, h);
       dxdt = U + Cg*k./sqrt(f^2 + Cg^2*dot(k, k, 2));
       dkdt = -([nablaU.u_x .* k(:,1) + nablaU.v_x .* k(:,2),nablaU.u_y .* k(:,1) + nablaU.v_y .* k(:,2)]);
       dydt = ode_xk2y(Npackets, dxdt, dkdt);
    end

    rayode = @odefun;
end

function dq=update(qk, K2, K_d2, beta, kx_, ky_)
   psik = -qk./(K_d2 + K2);
   psikx = 1i*kx_.*psik;
   psiky = 1i*ky_.*psik;
   qkx = 1i*kx_.*qk;
   qky = 1i*ky_.*qk;
   
   psix = k2g(psikx);
   psiy = k2g(psiky);
   qx = k2g(qkx);
   qy = k2g(qky);
   
   J = psix.*qy - psiy .* qx;
   %nu = 10;
   %hyperdiffusion = -nu*(1i*kx_).^2.*(1i*ky_).^2.*((1i*kx_).^6 + (1i*ky_).^6);
   dq = g2k(J) - beta*psikx;
end