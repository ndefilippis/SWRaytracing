test = zeros(1,frame);
for i=1:frame
    qk = g2k(q_save(:,:,i));
    psik = -qk./(K_d2 + K2);
    uk = -1i*ky_.*psik;
    vk =  1i*kx_.*psik;
    %KE = sum(sum(uk.*conj(uk) + vk.*conj(vk)));
    KEk = K2.*abs(psik).^2;
    KE = zeros(1,kmax);
    for j=1:kmax
       mask = ((j-1)^2 <= K2) .* (K2 < j^2);
       KE(j) = sum(mask(:).*KEk(:)); 
    end
    loglog(1:kmax, KE);
    hold on
    loglog(1:kmax, 1e-3*(1:kmax).^(-3), 'k--');
    hold off
    pause(1/30);
end
