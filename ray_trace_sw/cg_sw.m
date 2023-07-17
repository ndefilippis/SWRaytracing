function [C,omega,omega_abs,divC,gradomega] = cg_sw(k,l,C0,f,U,H)

% [C,omega,omega_abs,divC,gradomega] = cg_sw(k,l,C0,f,U,H)
%
% SW group velocity for wavenumber k, l, Coriolis f, and C0 =
% sqrt(g*H0). Optional H=H(:,:) is eta_g(x,y)/H0

% Returns C.x,C.y,omega, divergcence of C

%% ** unfixed half changes made!  Need local U, not spatially
%% dependent U, for absolute frequency.  You did this to check
%% whether absolute frequency is conserved


if nargin>5
    gH = C0^2*H;
else
    gH = C0^2;
end

%omega = sqrt( (f+vort/2)^2 + C0^2*(k^2+l^2));
omega = sqrt( f^2 + gH*(k^2+l^2) );
%omega_abs = omega + U.u*k + U.v*l;
omega_abs = abs(omega);
C.x = gH*k./omega;
C.y = gH*l./omega;

if nargout>2 && nargin>4
    divC = ( k*f*U.v - l*f*U.u - C.x.^2 - C.y.^2)./omega;
    gradomega.x = f*(k^2+l^2)*U.v./(2*omega);
    gradomega.y = -f*(k^2+l^2)*U.u./(2*omega);
end

