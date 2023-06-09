classdef RaytracingScheme
    methods(Abstract)
        psi = streamfunction(obj, x, y, t)
        u = U(obj, x, t)
        nablaU = grad_U(obj, x, t)
    end
    
    methods
        function nablaU_k = grad_U_times_k(obj, x, k, t)
            nablaU = obj.grad_U(x, t);            
            nablaU_k(:, 1) = nablaU.u_x .* k(:,1) + nablaU.v_x .* k(:,2);
            nablaU_k(:, 2) = nablaU.u_y .* k(:,1) + nablaU.v_y .* k(:,2);
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