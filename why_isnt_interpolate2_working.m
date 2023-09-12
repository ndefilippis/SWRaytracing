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
num_interp_points = 4;
%X_interp = linspace(-L/2 - num_interp_points*dx, L/2 + num_interp_points*dx, nx + 2*num_interp_points);
[XX, YY] = meshgrid(X);
XXi = [XX(:,end-num_interp_points:end-1) - L, XX, XX(:,2:num_interp_points+1) + L];
YYi = [YY(end-num_interp_points:end-1,:) - L; YY; YY(2:num_interp_points+1,:) + L];
XXi = [XXi(1:num_interp_points, :); XXi; XXi(end-num_interp_points+1:end, :)];
YYi = [YYi(:, 1:num_interp_points), YYi, YYi(:, end-num_interp_points+1:end)];
%[XXi, YYi] = meshgrid(X_interp);

%streamfunction = @(x, y, t) (0.5*(x.^2 + y.^2 - x.^4));
%scheme = DifferenceScheme(streamfunction);

q = read_field("analysis/pv", nx, nx, 1, [2000]);

scheme = SpectralScheme(L, nx, XXi, YYi, k2g(-g2k(q)./(K_d2 + K2)));

for i=1:100
t = linspace(0, L, nx);
slice_x = t;
slice_y = 0*t + L*i/100;

value1 = interpolate(slice_x, slice_y, scheme.U_field.u, dx, dx);
value2 = interpolate2(XXi, YYi, scheme.U_field.u, slice_x, slice_y, scheme.L); % HAD TO SWITCH X AND Y HERE

subplot(2, 1, 1);
contourf(XX, YY, scheme.U_field.u')
colorbar
hold on
plot(mod(slice_x, L), mod(slice_y, L), 'r');
subplot(2, 1, 2);
plot(t, value1, t, value2);
ylim([-0.5, 0.5]);
pause(1/15);
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