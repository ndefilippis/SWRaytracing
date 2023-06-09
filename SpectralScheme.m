classdef SpectralScheme < RaytracingScheme
    properties
        U_field, GradU_field, etag_field, L
    end
    methods
        function obj = SpectralScheme(Cg, f, L, nx, S)
            addpath C:/Users/ndefi/Downloads/rsw/
            
            obj.L = L;
            
            kmax = nx/2-1;
            [kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
            K2_ = kx_.^2 + ky_.^2;
            
            % In .mat file, variable Cg = sqrt(g*H0).  Here we'll
            % define gH0 = Cg^2, C0 = Cg
            gH0 = Cg^2;
            C0 = Cg;
            
            sig2_ = f^2 + gH0*K2_;
            
            uk   = g2k(S(:,:,1));
            vk   = g2k(S(:,:,2));
            etak = g2k(S(:,:,3));
            zetak  = 1i*(kx_.*vk - ky_.*uk);
            etagk = (f*etak-zetak).*f./sig2_;
            % zetagk = -gH0/f*etagk.*K2_;
            
            % Extract geostrophic velocities
            ugk = -1i*ky_.*(gH0/f*etagk);
            vgk = 1i*kx_.*(gH0/f*etagk);
            
            % Get gradients
            ugxk = 1i*kx_.*ugk;
            ugyk = 1i*ky_.*ugk;
            vgxk = 1i*kx_.*vgk;
            vgyk = 1i*ky_.*vgk;
            
            % Get fields in x-space
            obj.etag_field = k2g(etagk);
            obj.U_field.u = k2g(ugk);
            obj.U_field.v = k2g(vgk);
            
            obj.GradU_field.u_x = k2g(ugxk);
            obj.GradU_field.u_y = k2g(ugyk);
            obj.GradU_field.v_x = k2g(vgxk);
            obj.GradU_field.v_y = k2g(vgyk);
        end
        
        function psi = streamfunction(obj, x, y, t)
            dx = obj.L/size(obj.etag_field,1);
            xx = x(:);
            yy = y(:);
            psi = interpolate(xx,yy,obj.etag_field,dx,dx);
            psi = reshape(psi, size(x));
        end
        
        function u = U(obj, x, t)
            dx = obj.L/size(obj.U_field.u,1);
            xx = x(:, 1);
            yy = x(:, 2);
            u(:,1) = interpolate(xx,yy,obj.U_field.u,dx,dx);
            u(:,2) = interpolate(xx,yy,obj.U_field.v,dx,dx);
        end
        
        function nablaU = grad_U(obj, x, t)
            dx = obj.L/size(obj.GradU_field.u_x,1);
            xx = x(:, 1);
            yy = x(:, 2);
            nablaU.u_x = interpolate(xx,yy,obj.GradU_field.u_x,dx,dx);
            nablaU.u_y = interpolate(xx,yy,obj.GradU_field.u_y,dx,dx);
            nablaU.v_x = interpolate(xx,yy,obj.GradU_field.v_x,dx,dx);
            nablaU.v_y = interpolate(xx,yy,obj.GradU_field.v_y,dx,dx);
        end
    end
end
%
%
% function nablaU_k = grad_U_times_k(x, k, t)
% global GradU;
% dx = 1/size(GradU.u_x,1);
% xx = x(:, 1);
% yy = x(:, 2);
% U_x = interpolate(xx,yy,GradU.u_x,dx,dx);
% U_y = interpolate(xx,yy,GradU.u_y,dx,dx);
% V_x = interpolate(xx,yy,GradU.v_x,dx,dx);
% V_y = interpolate(xx,yy,GradU.v_y,dx,dx);
%
% nablaU_k(:, 1) = U_x .* k(:,1) + V_x .* k(:,2);
% nablaU_k(:, 2) = U_y .* k(:,1) + V_y .* k(:,2);
% end
%
% function eta=vorticity(x, t)
% nablaU = grad_U(x, t);
% eta = nablaU.v_x - nablaU.u_y;
% end
%
% function sigma=strain(x, t)
% nablaU = grad_U(x, t);
% sigma = sqrt((nablaU.u_x - nablaU.v_y).^2 + (nablaU.v_x + nablaU.u_y).^2);
% end
%
% function D = okuboWeiss(x, t)
% nablaU = grad_U(x, k, t);
% D = nablaU.v_y.^2 + nablaU.v_x.*nablaU.u_y;
% end
