function [x, k, t] =  ode_symplectic(x0, k0, dt, T, f, gH, scheme)
    Nsteps = floor(T / dt);
    x = zeros([Nsteps, size(x0, [2,3])]);
    k = zeros([Nsteps, size(k0, [2,3])]);
    t = zeros(Nsteps, 1);
    x(1, :, :) = x0;
    k(1, :, :) = k0;
    t(1) = 0;
    
    omega = @(k) (sqrt(f^2 + gH*dot(k, k, 2)));
    group_velocity = @(k) (gH * k./omega(k));
    
    function [x, k] = phi1(x0, k0, dt)
        x = x0 + dt * group_velocity(k0);
        k = k0;
    end

    function [x, k] = phi2(x0, k0, dt)
        x = x0 + dt * scheme.U(x0);
        k = k0 - dt * scheme.grad_U_times_k(x0, k0, 0);
    end
    
    for i=2:Nsteps
        [x0, k0] = leapfrog_method(x0, k0, dt, @phi1, @phi2);
        x(i, :, :) = x0;
        k(i, :, :) = k0;
        t(i) = (i-1)*dt;
    end
    
    
end

function [x, k] = leapfrog_method(x0, k0, dt, phi1, phi2)
        [x1, k1] = phi1(x0, k0, dt/2);
        [x2, k2] = phi2(x1, k1, dt);
        [x, k] = phi1(x2, k2, dt/2);
end

function [x, k] = yoshida(x0, k0, dt, phi1, phi2)
        two_root_three = nthroot(2, 3);
        w0 = two_root_three / (2 - two_root_three);
        w1 = 1/(2 - two_root_three);
        c1 = w1/2;
        c2 = (w0 + w1)/2;
        c3 = c2;
        c4 = c1;
        d1 = w1;
        d2 = w0;
        d3 = w1;
        [x1, k1] = phi1(x0, k0, dt/2);
        [x2, k2] = phi2(x1, k1, dt);
        [x, k] = phi1(x2, k2, dt/2);
end