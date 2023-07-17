x = linspace(-5, 5, 1000);
plt = plot(x, x);
xlim([-5, 5]);
ylim([0, 1]);
for i=1:1000
    t = i / 10;
    set(plt, 'YData', f(x, t));
    
    title(integrate(x, t));

    pause(1/60);
end

function y = f(x, t)
    L = 10;
    y = exp(-(mod(x + t + L/2, L) - L/2).^2);
end

function I=integrate(x, t)
    ff = @(x) f(x, t);
    I = integral(ff, -5, 5);
end