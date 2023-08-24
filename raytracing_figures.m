fmt = "yyyy-MM-dd HH_mm_ss";
filename = image_path + "QG_" + string(datetime("now"), fmt) + ".gif";

figure();
subplot(2, 1, 1);
p = scatter(x(:,1), x(:,2), 'k.');
axis image
hold on
p_solver = scatter(solver_x(1,:,1), solver_x(1,:,2), 'r.');
q = quiver(x(:,1), x(:,2), k(:,1), k(:,2), 0.5, 'k');
contour(XX, YY, scheme.streamfunction(XX, YY));
xlim([-L/2, L/2]);
ylim([-L/2, L/2]);
%legend("wavepacket")
xlabel("x");
ylabel("y");
grid on

subplot(2, 2, 3);
KK = scatter(k(:,1), k(:,2), 'k.');
xlabel("k");
ylabel("l");
grid on
axis image
xlim([-5, 5]);
ylim([-5, 5]);

subplot(2,2,4);
err_plot = plot(t_hist(1) * (f*Fr^2), error(1));
xlim([0, Tend * (f*Fr^2)]);
ylim([1.1*min(min(error)), 1.1*max(max(error))]);
%freq = sort(omega(squeeze(k_hist(1,:,:)), f, gH));
%W = plot(freq, energy(freq));

if(write_file)
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','DelayTime',0.05, 'Loopcount',inf);
end

anim = 0;
if(anim)
for i=1:40:Nsteps
       set(p, {'XData','YData'},{mod(x_hist(i,:,1) + L/2, L) - L/2, mod(x_hist(i,:,2) + L/2, L) - L/2});
       set(p_solver, {'XData', 'YData'},{mod(solver_x(i,:,1) + L/2, L) - L/2, mod(solver_x(i,:,2) + L/2, L) - L/2});
       set(q, {'XData', 'YData', 'UData', 'VData'}, {mod(x_hist(i,:,1) + L/2, L) - L/2, mod(x_hist(i,:,2) + L/2, L) - L/2, k_hist(i,:,1), k_hist(i,:,2)});
       set(KK, {'XData','YData'}, {k_hist(i,:,1), k_hist(i,:,2)});
       %freq = sort(omega(squeeze(k_hist(i,:,:)), f, gH));
       %set(W, {'XData', 'YData'}, {freq, energy(freq)});
       set(err_plot, {'XData','YData'}, {t_hist(1:i) * (f*Fr^2), error(1:i)})
       test = diff(error);
       
       %subplot(2, 2, 2);
       %subplot(1, 2, 2);
       %histogram(omega(k), 'Normalization', 'probability');
       %histogram(w.*sum(w' > w-5 & w' < w+5)');
       %scatter(omega(k)/f, energy(omega(k)));
       %xlim([0.1, 100]);
       %xlabel('\omega/f');
       %ylabel('e(\omega)');
       %set(gca,'xscale','log')
       %set(gca,'yscale','log')
       if(write_file)
           frame = getframe(1);
           im = frame2im(frame);
           [imind,cm] = rgb2ind(im,256);
           imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append');
       end
       pause(1/60);
end
end



figure()
title("Accumulation plot")
plot(k_hist(:,:,1), k_hist(:,:,2))
hold on
plot(k_hist(1,:,1), k_hist(1,:,2), 'k', 'LineWidth', 2.5)
hold off
xlabel("k_1");
ylabel("k_2");