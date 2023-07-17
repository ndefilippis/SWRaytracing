D = 40;
Npoints = 10000;
x = linspace(-D/2, D/2, Npoints);
dx = D / Npoints;

rng(10);
k = (5 * rand(3, 1) + 5) .* sign(rand(3, 1) - 0.5);
phase = D*rand(3, 1) - D/2;
f = 4;
g = 1;
H = 10;
gH = g * H;
omega = sqrt(f^2 + gH*k.^2);
Cg = gH*k./omega;

L = 2;
A = @(x,t) 0.3*exp(-(mod(x - Cg*t + phase + D/2, D) - D/2).^2/(2*L^2));

theta = @(x, t)  k*x - omega*t;
t = 0;
Nsteps = 2000;
dt = 0.1;
plt = plot(x, x, 'k');
hold on
top = plot(x, A(x, t), 'b');
bot = plot(x, -A(x, t), 'b');
hold off
ylim([-1, 1]);
xlim([-D/2, D/2]);
diff = zeros(Nsteps, 1);
for i=1:Nsteps
    eta_n = @(x) real(A(x, t).*exp(1i*theta(x, t)));
    eta = @(x) real(sum(A(x, t).*exp(1i*theta(x, t)), 1));
    u = @(x) real(sum(omega./(H*k).*A(x, t).*exp(1i*theta(x, t)), 1));
    v = @(x) real(sum(-1i*f./(H*k).*A(x, t).*exp(1i*theta(x, t)), 1));
    u_n = @(x) omega./(H*k).*eta_n(x);
    v_n = @(x) -1i*f./(H*k).*eta_n(x);
    energy_x = @(x) H*(u(x).*conj(u(x)) + v(x).*conj(v(x))) + g*(eta(x).*conj(eta(x)));
    energy_x_n = @(x) H*(u_n(x).*conj(u_n(x)) + v_n(x).*conj(v_n(x))) + g*(eta_n(x).*conj(eta_n(x)));
    energy_n = integral(energy_x_n, -D/2, D/2,'ArrayValued',true);
    energy = integral(energy_x, -D/2, D/2);
    t = t + dt;
    diff(i) = (energy - sum(energy_n)) / energy;
    set(plt, 'YData', eta(x));
    set(top, 'YData', sum(A(x, t), 1));
    set(bot, 'YData', -sum(A(x, t), 1));
    title(sprintf("%f\n%f", energy, sum(energy_n)));
    subtitle(energy_n);
    pause(1/60);
end

plot(diff);
