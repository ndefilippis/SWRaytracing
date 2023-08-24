function background_flow = grid_U(qk, K_d2, K2, kx_, ky_)
    psik = -qk./(K_d2 + K2);
    vk = 1i*kx_.*psik;
    uk = -1i*ky_.*psik;
    
    ukx = 1i*kx_.*uk;
    uky = 1i*ky_.*uk;
    vkx = 1i*kx_.*vk;
    vky = 1i*ky_.*vk;
    
    background_flow.u = k2g(uk);
    background_flow.v = k2g(vk);
    
    background_flow.ux = k2g(ukx);
    background_flow.uy = k2g(uky);
    background_flow.vx = k2g(vkx);
    background_flow.vy = k2g(vky);
end