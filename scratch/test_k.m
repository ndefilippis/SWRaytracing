H = 1;
g = 1;
f = 3;
k = linspace(0.1, 10);
L = 1;


t = 1;
x = 1;

omega = sqrt(f^2 + g*H*k.^2);
Cg = g*H*k/omega;

for t=1:100
    x = Cg*t + 0.1;
    u = (Cg*(x - Cg*t) - 1i*L^2*omega) ./ (H*(x - Cg*t) - 1i*k*H*L^2);
    v = -f*L^2 ./ (H*(x-Cg*t) - 1i*k*H*L^2);
    
    plot(k, real(u), 'k', k, imag(v), 'r');
    hold on
    plot(k, real(omega./(k*H)), 'k--', k, imag(-1i*f./(k*H)), 'r--');
    hold off;
    title("A = " + string(exp(-(x-Cg*t).^2/(2*L^2))));
    pause(1/20);
end

