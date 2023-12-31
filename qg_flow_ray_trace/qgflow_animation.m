nx = 256;
Npackets = 50;
f = 3;
Cg = 1;
shear_strength = 0.4;

L = 2*pi;
x = linspace(-L/2, L/2, nx);
[X, Y] = ndgrid(x, x);

kmax = nx/2-1;
[kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
K2 = kx_.^2 + ky_.^2;

K_d2 = f/Cg;

% Plotting options
t_background_save = read_field('data/job-36976465/run-2/pv_time');
frame = size(t_background_save, 2);
q_save = read_field('data/job-36976465/run-2/pv', nx, nx, 1, 1:frame);

t_packet_save = read_field('data/job-36976465/run-2/packet_time');
packet_frame = size(t_packet_save, 2);
x_save = read_field('data/job-36976465/run-2/packet_x', Npackets, 2, 1, 1:packet_frame);
k_save = read_field('data/job-36976465/run-2/packet_k', Npackets, 2, 1, 1:packet_frame);

T = max(t_packet_save);
packet_steps_per_mean_flow_step = round((t_background_save(2) - t_background_save(1)) / (t_packet_save(2) - t_packet_save(1)));
packet_frames = find(t_background_save >= t_packet_save(1));
if(length(packet_frames) == 0)
    packet_frame_start = 0;
else
    packet_frame_start = packet_frames(1);
end

do_write_video = true;
image_path = '../images/';
datefmt = "yyyy-MM-dd HH_mm_ss";
movie_filename = image_path + "QG_t_" + string(datetime("now"), datefmt) + ".avi";
v = VideoWriter(movie_filename);
q_max = max(abs(q_save), [], 'all');

% plotting variables for animation
pv_data = q_save(:,:,1);
packet_x_data = [];
packet_y_data = [];
packet_k_data = [];
packet_l_data = [];

% create plots
figure();
hold on
[~,qg_contour_plot] = contourf(X, Y, pv_data, 18,'LineColor','none');
qg_contour_plot.ZDataSource = 'pv_data';

packet_scatter_plot = scatter(packet_x_data, packet_y_data, 25, 'k.');
packet_scatter_plot.XDataSource = 'packet_x_data';
packet_scatter_plot.YDataSource = 'packet_y_data';
packet_quiver_plot = quiver(packet_x_data, packet_y_data, 0.3, 'k', 'LineWidth', 1);
packet_quiver_plot.XDataSource = 'packet_x_data';
packet_quiver_plot.YDataSource = 'packet_y_data';
packet_quiver_plot.UDataSource = 'packet_k_data';
packet_quiver_plot.VDataSource = 'packet_l_data';
hold off

axis equal
colormap(redblue)
clim([-q_max, q_max])
c = colorbar();
c.Label.String = "PV";
xlabel('X');
ylabel('Y');
xlim([-L/2, L/2]);
ylim([-L/2, L/2]);
title_text_array = {"";sprintf("\\fontsize{8}\\color{gray}\\rmf: %.1f, C_g: %.2f, T (1/f): %.3f", f, Cg, T*f); ""};
title(title_text_array);

open(v);
for i=1:frame
    background_flow = grid_U(g2k(q_save(:,:,i)), K_d2, K2, kx_, ky_, shear_strength);
    if(t_background_save(i) < t_packet_save(1) && mod(i, 3) ~= 1)
        continue
    end
    speed2 = background_flow.u.^2 + background_flow.v.^2;
    Umax = sqrt(max(speed2, [], 'all'));
    Umean = sqrt(mean(speed2, 'all'));
    if(Npackets > 0 && t_background_save(i) >= t_packet_save(1))
        for j = 1:2:packet_steps_per_mean_flow_step
            alpha = j / packet_steps_per_mean_flow_step;
            pv_data =  alpha * q_save(:,:,i) + (1 - alpha) * q_save(:,:,i - 1);
            save_index = (i - packet_frame_start) * packet_steps_per_mean_flow_step + j;
            if(save_index > packet_frame)
               break 
            end
            packet_x_data = squeeze(x_save(:,1,save_index));
            packet_y_data = squeeze(x_save(:,2,save_index));
            packet_k_data = squeeze(k_save(:,1,save_index));
            packet_l_data = squeeze(k_save(:,2,save_index));
            title_text_array{1} = sprintf("t = %6.3f (1/f)", t_packet_save(save_index)*f);
            title_text_array{3} = sprintf("\\fontsize{8}\\color{gray}\\rmU_{max}: %5.3f, Fr: %5.3f", Umax, Umean/Cg);
            title(title_text_array);
            refreshdata
            if(do_write_video)
                fig_frame = getframe(gcf);
                writeVideo(v,fig_frame);
            else
                pause(1/30);
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
        else
            pause(1/30);
        end
        
    end

end
close(v);