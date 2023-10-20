function background_flow = grid_U(qk, K_d2, K2, kx_, ky_, shear_strength)
    psik = -qk./(K_d2 + K2);
    vk = 1i*kx_.*psik;
    uk = -1i*ky_.*psik;
    
    ukx = 1i*kx_.*uk;
    uky = 1i*ky_.*uk;
    vkx = 1i*kx_.*vk;
    vky = 1i*ky_.*vk;
    
    background_flow.u = apply_3d(uk, @k2g) + shear_strength;
    background_flow.v = apply_3d(vk, @k2g);
    
    background_flow.ux = apply_3d(ukx, @k2g);
    background_flow.uy = apply_3d(uky, @k2g);
    background_flow.vx = apply_3d(vkx, @k2g);
    background_flow.vy = apply_3d(vky, @k2g);
end