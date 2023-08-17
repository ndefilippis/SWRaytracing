h = 0.2;
L = 2*pi;
dim = ceil(L/h);
x = linspace(0, L, dim);
y = linspace(0, L, dim);
[X, Y] = meshgrid(x, x);

dt = 0.3;
T = 100;
Nsteps = T/dt;

q = sin(X).*cos(Y);
q = q(:);
K_d2 = -1;

A = diagnostic_solver(K_d2, h, dim);
for step=1:Nsteps
    psi = A\q;
    
    psi_x = psi_x_deriv(psi, h, dim);
    psi_y = psi_y_deriv(psi, h, dim);
    
    B = prognostic_solver(psi_x, psi_y, h, dt, dim);
    q = B\q;
    contourf(reshape(q, [dim, dim]));
    
    pause(1/20);
    disp(step);
end


function A=diagnostic_solver(K_d2, h, dim)
    A = zeros(dim.^2);
    for i=1:dim
        for j=1:dim
            A(idx(i, j, dim), idx(i, j  , dim)) = K_d2 - 4/h.^2;
            A(idx(i, j, dim), idx(i-1, j, dim)) = 1/h.^2;
            A(idx(i, j, dim), idx(i+1, j, dim)) = 1/h.^2;
            A(idx(i, j, dim), idx(i, j-1, dim)) = 1/h.^2;
            A(idx(i, j, dim), idx(i, j+1, dim)) = 1/h.^2;
        end
    end
end

function B=prognostic_solver(psi_x, psi_y, h, dt, dim)
    B = zeros(dim.^2);
    for i=1:dim
        for j=1:dim
            B(idx(i, j, dim), idx(i, j, dim)) = 1;
            B(idx(i, j, dim), idx(i-1, j, dim)) = -psi_y(idx(i, j, dim))*dt/(2*h);
            B(idx(i, j, dim), idx(i+1, j, dim)) =  psi_y(idx(i, j, dim))*dt/(2*h);
            B(idx(i, j, dim), idx(i, j-1, dim)) =  psi_x(idx(i, j, dim))*dt/(2*h);
            B(idx(i, j, dim), idx(i, j+1, dim)) = -psi_x(idx(i, j, dim))*dt/(2*h);
        end
    end
end

function psi_x = psi_x_deriv(psi_solve, h, dim)
    psi = reshape(psi_solve, [dim,dim]);
    psi_x = (circshift(psi, 1, 2) - circshift(psi, -1, 2))/(2*h);
    psi_x = psi_x(:);
end

function psi_y = psi_y_deriv(psi_solve, h, dim)
    psi = reshape(psi_solve, [dim,dim]);
    psi_y = (circshift(psi, 1) - circshift(psi, -1))/(2*h);
    psi_y = psi_y(:);
end

% This should automatically implement a periodic domain
function n=idx(i, j, dim)
    n = ((mod(i-1, dim)+1)-1)*dim + mod(j-1, dim)+1;
end

