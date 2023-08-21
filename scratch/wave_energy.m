test = zeros(1,packet_frame);
for i=1:packet_frame
    omega = squeeze(sqrt(f^2 + Cg^2*dot(k_save(:,:,i), k_save(:,:,i), 2)));
    energy = omega.*sum(omega' > omega-3 & omega' < omega+3)';
    loglog(omega, energy);
    pause(1/30);
end
