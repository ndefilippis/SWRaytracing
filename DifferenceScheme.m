classdef DifferenceScheme < RaytracingScheme
    properties
        h, psi
    end
    methods
        
        % Expects a streamfunction of the form f(x, y, t);
        function obj = DifferenceScheme(stream)
            obj.h = nthroot(eps, 3);
            obj.psi = stream;
        end
        
        function psi = streamfunction(obj, x, y, t)
            if(nargin < 4)
                t = 0;
            end
            psi = obj.psi(x, y, t);
        end
        
        function u = U(obj, x, t)
            if(nargin < 3)
                t = 0;
            end
            xx = x(:, 1,:);
            yy = x(:, 2,:);
            u = zeros(size(x));
            u(:,2,:) = (obj.psi(xx + obj.h/2, yy, t) - obj.psi(xx - obj.h/2, yy, t))/obj.h;
            u(:,1,:) = -(obj.psi(xx, yy + obj.h/2, t) - obj.psi(xx, yy - obj.h/2, t))/obj.h;
        end
        
        function nablaU = grad_U(obj, x, t)
            if(nargin < 3)
                t = 0;
            end
            xx = squeeze(x(:, 1, :));
            yy = squeeze(x(:, 2, :));
            
            nablaU.v_x = (obj.psi(xx + obj.h, yy, t) - 2 * obj.psi(xx, yy, t) + obj.psi(xx - obj.h, yy, t))/obj.h/obj.h;
            nablaU.u_y = -(obj.psi(xx, yy + obj.h, t) - 2 * obj.psi(xx, yy, t) + obj.psi(xx, yy - obj.h, t))/obj.h/obj.h;
            
            nablaU.v_y = (obj.psi(xx + obj.h/2, yy + obj.h/2, t)...
                + obj.psi(xx - obj.h/2, yy - obj.h/2, t)...
                - obj.psi(xx - obj.h/2, yy + obj.h/2, t)...
                - obj.psi(xx + obj.h/2, yy - obj.h/2, t))/obj.h/obj.h;
            nablaU.u_x = -nablaU.v_y;
        end
    end
end

% function nablaU_k = grad_U_times_k(x, k, t)
% [psi_xx, psi_xy, psi_yy] = grad_U(x, t);
% nablaU_k(:, 1) = -psi_xy .* k(:,1) + psi_xx .* k(:,2);
% nablaU_k(:, 2) = -psi_yy .* k(:,1) + psi_xy .* k(:,2);
% end
%
% function eta=vorticity(x, t)
% [psi_xx, ~, psi_yy] = grad_U(x, t);
% eta = psi_xx + psi_yy;
% end
%
% function sigma=strain(x, t)
% [psi_xx, psi_xy, psi_yy] = grad_U(x, t);
% sigma = 2*sqrt(psi_xy.^2 + (psi_xx - psi_yy).^2);
% end
%
% function D = okuboWeiss(x, t)
% [psi_xx, psi_xy, psi_yy] = grad_U(x, t);
% D = psi_xy.^2 - psi_xx.*psi_yy;
% end