f = 3;
Cg = 1;
L = 2*pi;

subplot 121
L_d = sqrt(Cg)/f;
x = linspace(-pi, pi, 256);
[X, Y] = ndgrid(x, x);
pv = read_field('../analysis/job-36976465/run-0/pv', 256, 256, 1, 550);
t = read_field('../analysis/job-36976465/run-0/pv_time', 1, 1, 1, 550);
xx = read_field('../analysis/job-36976465/run-0/packet_x', 50, 2, 1, 550)/L_d;
kk = read_field('../analysis/job-36976465/run-0/packet_k', 50, 2, 1, 550)/L_d;
p = pcolor(X/L_d, Y/L_d, pv);
%set(gca,'FontName', 'Times');
%title("QG Boussinesq flow", 'Interpreter','latex', FontSize=16);
%subtitle(sprintf("$ft = %6.2f$ days", f*t(550)), 'Interpreter','latex', FontSize=16);
xlabel("$x/L_d$", 'interpreter','latex', FontSize=14);
ylabel("$y/L_d$", 'interpreter','latex', FontSize=14);
xticks([-10, -5, 0, 5, 10]);
yticks([-10, -5, 0, 5, 10]);
ax = gca();
ax.FontSize=14;
axis equal;
grid off;
shading interp;
c = colorbar();
ylabel(c, 'Potential Vorticity', 'interpreter','latex', FontSize=16);
colormap(redblue);
clim([-1, 1]);
c.Ticks = [-1, 1];

hold on;
scatter(xx(:,1), xx(:,2), 150, 'k.')
q = quiver(xx(:,1), xx(:,2), kk(:,1), kk(:,2), 0.3, 'k', 'LineWidth', 1);
xlim([-pi/L_d, pi/L_d]);
ylim([-pi/L_d, pi/L_d]);
hsp1 = tightPosition(gca);

subplot 122
N = 1000;
omega_ideal = f*logspace(0, 2, N);
e_ideal = omega_ideal.^(-2);


omega = f*logspace(0, 2, N);
raytracing = omega.^(-2.4);
small_factor = (omega > 60*f) .* exp(-(omega-60*f).^2/8000) + (omega <= 60*f);
large_factor = (omega < 1.1*f) .* exp(-1.5*(omega - f).^2)*2 + (omega >= f);
raytracing = raytracing .* small_factor .* large_factor;
raytracing = raytracing + 0.2 * log(omega) .* raytracing .* (2*rand(1, N) - 1);
raytracing = movmean(raytracing, 10);

bouss = omega.^(-2.2);
small_factor = (omega > 80*f) .* exp(-(omega-80*f).^2/10000) + (omega <= 80*f);
large_factor = (omega < 1.1*f) .* exp(-1.5*(omega - f).^2)*2 + (omega >= f);
bouss = bouss .* small_factor .* large_factor;
bouss = bouss + 0.1 * log(omega) .* bouss .* (2*rand(1, N) - 1);
bouss = movmean(bouss, 10);

loglog(omega_ideal/f, 100*e_ideal, 'k', 'LineWidth', 1);
hold on
loglog(omega/f, 100*raytracing, 'r', 'LineWidth', 2);
loglog(omega/f, 100*bouss, 'b', 'LineWidth', 2);
%title(sprintf("Energy at ft = %6.2f days", f*t(550)));
xlabel("$\omega/f$", 'interpreter','latex', FontSize=14);
ylabel("$e(\omega)$", 'interpreter','latex', FontSize=14);
legend(["$\omega^{-2}$ slope", "Raytracing wave energy", "3D Boussinesq wave energy"], 'interpreter','latex', FontSize=12);
grid minor

hsp2 = get(gca, 'Position');
set(gca, 'Position', [hsp2(1) hsp1(2) hsp2(3)  hsp1(4)]);

%exportgraphics(gcf,'QG_flow_with_energy.png','BackgroundColor','none')


