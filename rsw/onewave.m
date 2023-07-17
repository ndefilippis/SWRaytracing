function [u,v,eta] = onewave(eta0,f,C0,k,l,omega,theta);

% theta = k*x_ + l*y_ + phi - (omega + kU + lV)*t
% omega = omegasign*sqrt(f^2+C0^2*K2);

eta = eta0*cos(theta);
u = eta0*(k*omega.*cos(theta)-l*f*sin(theta))/(k^2+l^2);
v = eta0*(l*omega.*cos(theta)+k*f*sin(theta))/(k^2+l^2);
