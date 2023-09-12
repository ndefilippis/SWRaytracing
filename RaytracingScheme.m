classdef RaytracingScheme
    methods(Abstract)
        psi = streamfunction(obj, x, y, t)
        u = U(obj, x, t)
        nablaU = grad_U(obj, x, t)
    end
    
    methods
        function nablaU_k = grad_U_times_k(obj, x, k, t)
            nablaU = obj.grad_U(x, t);        
            kk = k(:,1,:);
            ll = k(:,2,:);
            nablaU_k = zeros(size(k));
            nablaU_k(:, 1,:) = reshape(nablaU.u_x .* kk(:) + nablaU.v_x .* ll(:), size(kk));
            nablaU_k(:, 2,:) = reshape(nablaU.u_y .* kk(:) + nablaU.v_y .* ll(:), size(ll));;
        end
        
        function eta=vorticity(obj, x, t)
            nablaU = obj.grad_U(x, t);
            eta = nablaU.v_x - nablaU.u_y;
        end
        
        function sigma=strain(obj, x, t)
            nablaU = obj.grad_U(x, t);
            sigma = sqrt((nablaU.u_x - nablaU.v_y).^2 + (nablaU.v_x + nablaU.u_y).^2);
        end
        
        function D = okuboWeiss(obj, x, t)
            nablaU = obj.grad_U(x, k, t);
            D = nablaU.v_y.^2 + nablaU.v_x.*nablaU.u_y;
        end
    end
end