subplot(2, 1, 1);

t = linspace(0, 2*pi);
k_0 = 3;
k_vector = k_0*[cos(t)', sin(t)'];
U_dot_k = U(:,1)*k_vector(:,1)' + U(:,2)*k_vector(:,2)';
U_dot_k = U_dot_k(:);

omega_0 = sqrt(f^2 + Cg^2*k_0^2);
omega_abs = omega_0 + U_dot_k;
histogram(omega_abs, 'Normalization', 'pdf');
xlim([2, 6]);
ylabel("pdf");
title("Theoretical distribution of \omega");

subplot(2, 1, 2);

gap = 1;
samples = w(1:gap:end,:);
samples = samples(:);
histogram(samples, 'Normalization', 'pdf');
xlim([2, 6]);
ylabel("pdf");
title("Experimental distribution of \omega");
xlabel("\omega");