addpath ../qg_flow_ray_trace
num_runs = 9;
f_array = [];
Cg_array = [];
Ug_array = [];
legend_objs = [];
figure();
param_index = 1;
job_number = "job-37011720";
for i=16:16
    directory = job_number + "/run-" + num2str(i) + "/";
    
    [fid,status]=fopen(directory + "run.log");
    if(status)
        disp(status);
        continue;
    end
    resolution_cell = textscan(fid,'Resolution: %dx%d',1,'delimiter','\n', 'headerlines', 10);    
    Npackets_cell = textscan(fid, "Number of packets: %d",1);
    f_cell = textscan(fid, "Coriolis parameter: %f",'delimiter','\n', 'headerlines', 7);
    Cg_cell = textscan(fid,  "Group velocity: %f", 1);
    Ug_cell = textscan(fid, "Background velocity (parameter,computed): (%f,%f)", 1);
    resolution = resolution_cell{1};
    Npackets = Npackets_cell{1};
    f = f_cell{1};
    Cg = Cg_cell{1};
    Ug = Ug_cell{1};
    
    t_packet_save = read_field(directory + "packet_time");
    packet_frame = size(t_packet_save, 2);
    x_save = read_field(directory + "packet_x", Npackets, 2, 1, 1:packet_frame);
    k_save = read_field(directory + "packet_k", Npackets, 2, 1, 1:packet_frame);
    omega = squeeze(sqrt(f^2 + Cg^2.*dot(k_save, k_save, 2)));
    
    offset = 500;
    times = [1, 1000, 30000, length(t_packet_save)-offset];
    style = ['k', 'r', 'b', 'm'];
    bins = 300;
    edges = linspace(0, max(omega, [], "all"), bins);
    center = (edges(2:end) + edges(1:end-1))/2;
    
    scatter(omega(1,1)/f, omega(1, 1)*length(omega(:,1))*2*offset, 500, 'k.');
    hold on
    for t_index=2:length(times)
        w = omega(:,times(t_index)-offset:1:times(t_index)+offset);
        w = sort(w(:));
        distribution = histcounts(w, edges);
        %energy = w .* sum((w'-0.5 < w) .* (w < w'+0.5))';
        energy = center .* distribution;
        p = loglog(center/f, energy, style(t_index), "LineWidth", 2);
        %p = loglog(w, energy, style(t_index), "LineWidth", 2);
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    
    
    %if(f(1) ~= 1) 
    %    continue;  
    %end
    
    [color, line_style] = assign_style(f, Cg, Ug);

    mean_omega = mean(omega) / f;
    p = plot(t_packet_save * f, mean_omega);
    loglog(w/f, energy, color + line_style);
    if(true)
        Cg_array(param_index) = Cg;
        f_array(param_index) = f;
        Ug_array(param_index) = Ug;
        legend_objs(param_index) = p;
        param_index = param_index + 1;
    end
    hold on
    %p = plot(t_packet_save * f, mean_omega);
    
    title_string = "Energy versus omega (f=" + f + ",Cg=" + Cg + ",Ug=" + Ug + ",\omega_0/f=" + round(omega(1,1)/f) + ")";
    title(title_string);
    xlabel("\omega/f")
    ylabel("e(\omega)");
    %title("Energy distributions (f = " + f(1) + ")");
    %xlabel("\omega (f)");
    %ylabel("e(\omega)");
    %legend(legend_objs, "Ug = " + Ug_array + ", Cg = " + Cg_array + ", Fr = " + Ug_array ./ Cg_array);
    legend("ft = " + ((t_packet_save(times) - t_packet_save(1))*f));
    hold off
    
    filename = job_number + "_" + replace(replace(title_string, "\omega_0/f", "w0"), " ", "_") + ".png";
    saveas(gcf,filename)
end


function [color, line_style]=assign_style(f, Cg, Ug)
    color = "k";
    line_style = "-";
    if(Cg == 0.5)
        color = "r";
    elseif(Cg == 1.0)
        color = "b";
    end
    if(Ug == 0.2)
        line_style = "--";
    end
end