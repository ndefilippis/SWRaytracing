L = 1;
Cg = 0;
cor_f = 1;
gH = 1;
A = @(x, t) exp(-(x-Cg*t).^2/(2*L^2));
%A = @(x, t) (abs(x) < L)
t = 0;

ks = logspace(0, 1, 1000);
value = zeros(size(ks));
value2 = zeros(size(ks));
for i=1:length(ks)
    k = ks(i);
    omega = sqrt(cor_f^2 + gH*k.^2);
    f = @(x) (A(x, t).*real(exp(1i*(k*x - omega*t)))).^2;
    value(i) = integral(f, -inf, inf);
end
plot(ks, value, 'k');
hold on
theoretical_E = sqrt(pi)/2.*(exp(-ks.^2*L^2) + 1)*L;
plot(ks, theoretical_E, 'b--');

