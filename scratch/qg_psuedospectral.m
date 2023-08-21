addpath ../rsw/

% Set up domain
L = 2*pi;
nx = 128;
h = L/nx;
x = linspace(-L/2, L/2, nx);
y = linspace(-L/2, L/2, nx);
[X, Y] = meshgrid(x, x);

kmax = nx/2-1;
[kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
K2 = kx_.^2 + ky_.^2;

% Simulation parameters
beta = 0;
f = 3;
Cg = 1;
K_d2 = f/Cg;

% Set up wave packets
Npackets = 30;

packet_x = zeros(Npackets, 2);
packet_k = zeros(Npackets, 2);

for i=1:Npackets
   packet_k(i, :) = 1.1 * f * [cos(2*pi*i/Npackets), sin(2*pi*i/Npackets)]; 
   packet_x(i, :) = L*(rand(1, 2)) - L/2;
end
packet_steps_per_simulation_step = 10;

a_g = 0.1;
q = initial_q(X, Y, a_g, K_d2);
qk = g2k(q);

% Find CFL condition based on inital PV
% Also compute Froude number from inital flow
flow = grid_U(qk, K_d2, K2, kx_, ky_);
speed2 = flow.u.^2 + flow.v.^2;
U0 = sqrt(max(speed2(:)));
Fr = U0/Cg;
fprintf("Froude Number: %f\n", Fr);

dt = 0.3*h/U0;
packet_dt = dt/packet_steps_per_simulation_step;
T = 1000/f;
Nsteps = ceil(T/dt);
steps_per_save = 20; % Has to be bigger than 3 or it wont save the initial AB steps

% Npacket_steps = Nsteps * packet_steps_per_simulation_step;
packet_steps_per_save = 1;
packet_delay = 0.75*T;
packet_step_start = ceil(packet_delay / T * Nsteps);
packet_frame_start = ceil(packet_delay / T * Nsteps / steps_per_save);

q_save = zeros(nx, nx, floor(Nsteps / steps_per_save));
x_save = zeros(Npackets, 2, floor((Nsteps - packet_step_start) / packet_steps_per_save));
k_save = zeros(Npackets, 2, floor((Nsteps - packet_step_start) / packet_steps_per_save));
t_background_save = zeros(1, floor(Nsteps / steps_per_save));
t_packet_save = zeros(1, floor((Nsteps - packet_step_start) / packet_steps_per_save));

frame = 1;
packet_frame = 1;
q_save(:, :, frame) = q;
x_save(:, :, packet_frame) = packet_x;
k_save(:, :, packet_frame) = packet_k;

AB_order = 3;
Q_n_minus = zeros(2*kmax+1, kmax+1, AB_order - 1);
% Initialize Adams-Bashforth inital steps using:
% 1. Forward Euler
% 2. Second order AB
Q_n_minus(:,:,2) = update(qk, K2, K_d2, beta, kx_, ky_);
qk = qk + dt*Q_n_minus(:,:,2);

Q_n_minus(:,:,1) = update(qk, K2, K_d2, beta, kx_, ky_);
qk = qk + dt/2*(3*Q_n_minus(:,:,1) - Q_n_minus(:,:,2));

tic
fprintf("Simulation progress:  0.00%%")
Ef = filter(kx_, ky_, h);
t = 2*dt;
S = sparsity(Npackets);
for step=4:Nsteps
   prev_qk = qk;
   [dq, Qn] = AB3(qk, K2, K_d2, beta, kx_, ky_, dt, Q_n_minus);
   t = t + dt;
   Q_n_minus(:,:,2) = Q_n_minus(:,:,1);
   Q_n_minus(:,:,1) = Qn;
   qk = qk + dq;
   qk = Ef  .* qk;% At some point, try to get hyperdiffusion working here
   
   % Do wavepacket advection
   if(t > packet_delay)
       background_flow1 = grid_U(prev_qk, K_d2, K2, kx_, ky_);
       background_flow2 = grid_U(qk, K_d2, K2, kx_, ky_);
       ray_ode = generate_raytracing_ode(background_flow1, background_flow2, Npackets, f, Cg, dt, h);
       y0 = ode_xk2y(Npackets, packet_x, packet_k);
       opts = odeset('Jpattern', S);
       [t_steps, solver_y] = ode15s(ray_ode, [0, dt], y0, opts);
       [packet_x, packet_k] = ode_y2xk(Npackets, solver_y(end,:)');
   end
   
   if(t > packet_delay && mod(step, packet_steps_per_save) == 0)
      packet_frame = packet_frame + 1; 
      x_save(:, :, packet_frame) = mod(packet_x + L/2, L) - L/2;
      k_save(:, :, packet_frame) = packet_k;
      t_packet_save(packet_frame) = t;
   end
   
   if(mod(step, steps_per_save) == 0)
       frame = frame + 1;
       q_save(:,:,frame) = k2g(qk);
       t_background_save(frame) = t;
   end
   if mod(step, 51) == 0
        fprintf("\b\b\b\b\b\b\b% 6.2f%%", step/Nsteps*100)
   end
end
fprintf("\b\b\b\b\b\b\b100.00%%\n");
time_elapsed = toc;

% Plotting options

do_write_video = true;
image_path = '../images/';
datefmt = "yyyy-MM-dd HH_mm_ss";
movie_filename = image_path + "QG_t_" + string(datetime("now"), datefmt);
v = VideoWriter(movie_filename);
q_max = max(abs(q_save), [], 'all');

% plotting variables for animation
pv_data = q_save(:,:,1);
packet_x_data = [];
packet_y_data = [];

% create plots
figure();
hold on
[~,qg_contour_plot] = contourf(X, Y, pv_data, 18,'LineColor','none');
qg_contour_plot.ZDataSource = 'pv_data';

packet_scatter_plot = scatter(packet_x_data, packet_y_data, 25, 'k.');
packet_scatter_plot.XDataSource = 'packet_x_data';
packet_scatter_plot.YDataSource = 'packet_y_data';
hold off

axis image
colormap(redblue)
caxis([-q_max, q_max])
c = colorbar();
c.Label.String = "PV";
xlabel('X');
ylabel('Y');
title_text_array = {"";sprintf("\\fontsize{8}\\color{gray}\\rmf: %.1f, C_g: %.2f, T (1/f): %.3f, time elapsed: %.2f sec", f, Cg, T*f, time_elapsed); ""};
title(title_text_array);

open(v);
for i=1:frame
    background_flow = grid_U(g2k(q_save(:,:,i)), K_d2, K2, kx_, ky_);
    speed2 = background_flow.u.^2 + background_flow.v.^2;
    Umax = max(speed2, [], 'all');
    if(i > packet_frame_start)
        for j = 1:steps_per_save / packet_steps_per_save
            alpha = packet_steps_per_save * j / steps_per_save;
            pv_data =  alpha * q_save(:,:,i) + (1 - alpha) * q_save(:,:,i - 1);
            save_index = (i - packet_frame_start - 1) * steps_per_save / packet_steps_per_save + j;
            if(save_index > packet_frame)
               break 
            end
            packet_x_data = squeeze(x_save(:,1,save_index));
            packet_y_data = squeeze(x_save(:,2,save_index));
            title_text_array{1} = sprintf("t = %6.3f (1/f)", t_packet_save(save_index)*f);
            title_text_array{3} = sprintf("\\fontsize{8}\\color{gray}\\rmU_{max}: %5.3f, Fr: %5.3f", Umax, Umax/Cg);
            title(title_text_array);
            refreshdata
            if(do_write_video)
                fig_frame = getframe(gcf);
                writeVideo(v,fig_frame);
            end
        end
    else
        pv_data =  q_save(:,:,i);
        title_text_array{1} = sprintf("t = %6.3f (1/f)", t_background_save(i)*f);
        title_text_array{3} = sprintf("\\fontsize{8}\\color{gray}\\rmU_{max}: %5.3f, Fr: %5.3f", Umax, Umax/Cg);
        title(title_text_array);
        refreshdata
        if(do_write_video)
            fig_frame = getframe(gcf);
            writeVideo(v,fig_frame);
        end
    end

end
close(v);

function q=initial_q(X, Y, a_g, K_d2)
    % Inital background PV
    % set as a ring of intermediate wavenumbers
    q = 0*X;
    phase = 2*pi*rand(31, 31);
    total_eta = 0;
    for k=-13:13
        for l=-13:13
            if 5^2 < k^2 + l^2 <= 8^2
                % Keep track of total eta to normalize after the face
                eta = real(exp(1i*(k*X + l*Y + phase(k+16, l+16))));
                total_eta = total_eta + eta;

                % The constant here is to satisfy geostrophy
                q = q + 1/K_d2*(k^2+l^2)*eta;
            end
        end
    end
    q = a_g * 1/max(max(total_eta)) * q;
end

function Ef=filter(kx_, ky_, h)
    Ef = ones(size(kx_));
    kstar = sqrt((kx_*h).^2 + (ky_*h).^2);
    
    kc = 0.75*pi;
    const = log(1e-15)/(0.25*pi)^4;
    result = exp(const * (kstar - kc).^4);
    Ef(kstar >= kc) = result(kstar >= kc);
end

function update=rk4(dt, k, x, f)
    k1 = f(k, x);
    k2 = f(k + dt * k1/2, x + dt * k1/2);
    k3 = f(k + dt * k2/2, x + dt * k2/2);
    k4 = f(k + dt * k3, x + dt * k3);
    
    update = dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

%function dq = RK4(qk, K2, K_d2, beta, kx_, ky_, dt)
%   k1 = update(qk, K2, K_d2, beta, kx_, ky_);
%   k2 = update(qk + k1*dt/2, K2, K_d2, beta, kx_, ky_);
%   k3 = update(qk + k2*dt/2, K2, K_d2, beta, kx_, ky_);
%   k4 = update(qk + k3*dt, K2, K_d2, beta, kx_, ky_);
%   dq = dt/6 * (k1 + 2*k2 + 2*k3 + k4);
%end

function [dq, Qn] = AB3(qk, K2, K_d2, beta, kx_, ky_, dt, Q_n_minus)
    Qn = update(qk, K2, K_d2, beta, kx_, ky_);
    dq = dt/12*(23*Qn - 16*Q_n_minus(:,:,1) + 5*Q_n_minus(:,:,2));
end

function background_flow = grid_U(qk, K_d2, K2, kx_, ky_)
    psik = -qk./(K_d2 + K2);
    vk = 1i*kx_.*psik;
    uk = -1i*ky_.*psik;
    
    ukx = 1i*kx_.*uk;
    uky = 1i*ky_.*uk;
    vkx = 1i*kx_.*vk;
    vky = 1i*ky_.*vk;
    
    background_flow.u = k2g(uk);
    background_flow.v = k2g(vk);
    
    background_flow.ux = k2g(ukx);
    background_flow.uy = k2g(uky);
    background_flow.vx = k2g(vkx);
    background_flow.vy = k2g(vky);
  
end

function [U, nablaU] = interpolate_U(background_flow1, background_flow2, alpha, x, h)
    xx = x(:, 1);
    yy = x(:, 2);
    
    U1(:,1) = interpolate(xx,yy,background_flow1.u,h,h);
    U1(:,2) = interpolate(xx,yy,background_flow1.v,h,h);
    U2(:,1) = interpolate(xx,yy,background_flow2.u,h,h);
    U2(:,2) = interpolate(xx,yy,background_flow2.v,h,h);
    
    nablaU1.u_x = interpolate(xx,yy,background_flow1.ux,h,h);
    nablaU1.u_y = interpolate(xx,yy,background_flow1.uy,h,h);
    nablaU1.v_x = interpolate(xx,yy,background_flow1.vx,h,h);
    nablaU1.v_y = interpolate(xx,yy,background_flow1.vy,h,h);
    nablaU2.u_x = interpolate(xx,yy,background_flow2.ux,h,h);
    nablaU2.u_y = interpolate(xx,yy,background_flow2.uy,h,h);
    nablaU2.v_x = interpolate(xx,yy,background_flow2.vx,h,h);
    nablaU2.v_y = interpolate(xx,yy,background_flow2.vy,h,h);
    
    U = (1 - alpha) * U1 + alpha * U2;
    nablaU.u_x = (1 - alpha) * nablaU1.u_x + alpha * nablaU2.u_x;
    nablaU.u_y = (1 - alpha) * nablaU1.u_y + alpha * nablaU2.u_y;
    nablaU.v_x = (1 - alpha) * nablaU1.v_x + alpha * nablaU2.v_x;
    nablaU.v_y = (1 - alpha) * nablaU1.v_y + alpha * nablaU2.v_y;
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

function result=update(qk, K2, K_d2, beta, kx_, ky_)
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
   result = g2k(J) - beta*psikx;
end