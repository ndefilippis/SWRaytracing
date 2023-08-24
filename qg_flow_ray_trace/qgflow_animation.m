nx = 128;
Npackets = 40;
f = 3;
Cg = 1;

L = 2*pi;
x = linspace(-L/2, L/2, nx);
[X, Y] = meshgrid(x, x);

kmax = nx/2-1;
[kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
K2 = kx_.^2 + ky_.^2;

K_d2 = f/Cg;

% Plotting options
t_background_save = read_field('pv_time');
frame = size(t_background_save, 2);
q_save = read_field('pv', nx, nx, 1, 1:frame);

t_packet_save = read_field('packet_time');
packet_frame = size(t_packet_save, 2);
x_save = read_field('packet_x', Npackets, 2, 1, 1:packet_frame);
k_save = read_field('packet_k', Npackets, 2, 1, 1:packet_frame);

T = max(t_packet_save);
packet_steps_per_mean_flow_step = round((t_background_save(2) - t_background_save(1)) / (t_packet_save(2) - t_packet_save(1)));
packet_frames = find(t_background_save >= t_packet_save(1));
packet_frame_start = packet_frames(1);

do_write_video = false;
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
title_text_array = {"";sprintf("\\fontsize{8}\\color{gray}\\rmf: %.1f, C_g: %.2f, T (1/f): %.3f", f, Cg, T*f); ""};
title(title_text_array);

open(v);
for i=1:frame
    background_flow = grid_U(g2k(q_save(:,:,i)), K_d2, K2, kx_, ky_);
    speed2 = background_flow.u.^2 + background_flow.v.^2;
    Umax = max(speed2, [], 'all');
    if(t_background_save(i) >= t_packet_save(1))
        for j = 1:packet_steps_per_mean_flow_step
            alpha = j / packet_steps_per_mean_flow_step;
            pv_data =  alpha * q_save(:,:,i) + (1 - alpha) * q_save(:,:,i - 1);
            save_index = (i - packet_frame_start) * packet_steps_per_mean_flow_step + j;
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