addpath ../qg_flow_ray_trace
num_runs = 9;
f_array = [];
Cg_array = [];
Ug_array = [];
legend_objs = [];
figure();
param_index = 1;
for i=1:num_runs
    directory = "job-36976465/run-" + num2str(i-1) + "/";
    
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
    w = sort(omega(:,end));
    energy = w .* sum((w'-1 < w) .* (w < w'+1))';
    %if(f(1) ~= 1) 
    %    continue;  
    %end
    
    [color, line_style] = assign_style(f, Cg, Ug);

    mean_omega = mean(omega) / f;
    p = plot(t_packet_save * f, mean_omega);
    %loglog(w/f, energy, color + line_style);
    if(true)
        Cg_array(param_index) = Cg;
        f_array(param_index) = f;
        Ug_array(param_index) = Ug;
        legend_objs(param_index) = p;
        param_index = param_index + 1;
    end
    hold on
end
title("Packet-averaged wavenumber versus time (f = " + f_array(1) + ")");
xlabel("t (1/f)")
ylabel("\omega (f)");
%title("Energy distributions (f = " + f(1) + ")");
%xlabel("\omega (f)");
%ylabel("e(\omega)");
legend(legend_objs, "Ug = " + Ug_array + ", Cg = " + Cg_array + ", Fr = " + Ug_array ./ Cg_array);

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