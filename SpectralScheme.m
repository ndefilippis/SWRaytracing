classdef SpectralScheme < RaytracingScheme
    properties
        U_field, GradU_field, psi_field, X_interp, Y_interp, L
    end
    methods
        function obj = SpectralScheme(L, nx, psi_field)
            addpath ./rsw/
            addpath ./ray_trace_sw/
            
            obj.L = L;
            
            kmax = nx/2-1;
            [kx_,ky_] = ndgrid(-kmax:kmax,0:kmax);
            K2_ = kx_.^2 + ky_.^2;
            
            psik   = g2k(psi_field);            
            % Extract geostrophic velocities
            ugk = -1i*ky_.*psik;
            vgk = 1i*kx_.*psik;
            
            % Get gradients
            ugxk = 1i*kx_.*ugk;
            ugyk = 1i*ky_.*ugk;
            vgxk = 1i*kx_.*vgk;
            vgyk = 1i*ky_.*vgk;
            
            % Get fields in x-space
            obj.psi_field = k2g(psik);
            obj.U_field.u = k2g(ugk);
            obj.U_field.v = k2g(vgk);
            
            obj.GradU_field.u_x = k2g(ugxk);
            obj.GradU_field.u_y = k2g(ugyk);
            obj.GradU_field.v_x = k2g(vgxk);
            obj.GradU_field.v_y = k2g(vgyk);
        end
        
        function psi = streamfunction(obj, x, y, t)
            dx = obj.L/size(obj.psi_field, 1);
            psi = interpolate(x,y,obj.psi_field,dx,dx);
            %psi = interpolate2(obj.X_interp,obj.Y_interp,obj.psi_field,x,y,obj.L);
            psi = reshape(psi, size(x));
        end
        
        function u = U(obj, x, t)
            dx = obj.L/size(obj.psi_field, 1);
            xx = x(:,1,:);
            yy = x(:,2,:);
            u = zeros(size(x));
            u(:,1,:) = reshape(interpolate(xx(:),yy(:),obj.U_field.u,dx,dx), size(xx));
            u(:,2,:) = reshape(interpolate(xx(:),yy(:),obj.U_field.v,dx,dx), size(xx));
            %u(:,1,:) = interpolate2(obj.X_interp,obj.Y_interp,obj.U_field.u,xx,yy,obj.L);
            %u(:,2,:) = interpolate2(obj.X_interp,obj.Y_interp,obj.U_field.v,xx,yy,obj.L);
        end
        
        function nablaU = grad_U(obj, x, t)
            dx = obj.L/size(obj.psi_field, 1);
            xx = x(:,1,:);
            yy = x(:,2,:);
            nablaU.u_x = interpolate(xx(:),yy(:),obj.GradU_field.u_x,dx,dx);
            nablaU.u_y = interpolate(xx(:),yy(:),obj.GradU_field.u_y,dx,dx);
            nablaU.v_x = interpolate(xx(:),yy(:),obj.GradU_field.v_x,dx,dx);
            nablaU.v_y = interpolate(xx(:),yy(:),obj.GradU_field.v_y,dx,dx);
            %nablaU.u_x = interpolate2(obj.X_interp,obj.Y_interp,obj.GradU_field.u_x,xx,yy,obj.L);
            %nablaU.u_y = interpolate2(obj.X_interp,obj.Y_interp,obj.GradU_field.u_y,xx,yy,obj.L);
            %nablaU.v_x = interpolate2(obj.X_interp,obj.Y_interp,obj.GradU_field.v_x,xx,yy,obj.L);
            %nablaU.v_y = interpolate2(obj.X_interp,obj.Y_interp,obj.GradU_field.v_y,xx,yy,obj.L);
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
