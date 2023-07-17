function [xp] = advect1d(xpin,u,x,dt)
    
% [xp] = advect1d(xpin,u,x,dt)
% particle position vector for particle m: X_m(t)
% m = 1:np
    
% advect with RK4

nx = length(u);
L = 2*x(end)-x(end-1)-x(1);  % 2*pi
dx = L/nx;

% since u is periodic on x, extend it a few points in either
% direction to ensure we're never interpolating at end points

u = u(:)';  % make sure its a row vector
x = x(:)';

ne = 10;
ub = [u(end-ne+1:end) u u(1:ne)];
xb = [x(1)-dx*[ne:-1:1] x x(end)+dx*[1:ne]];

% Get position of each particle on [x(1) x(end)]
xp0 = mod(xpin-x(1),L)+x(1);  

% ... and also number of L-lengths from [0,L)
xoff = L*floor((xpin-x(1))/L);

xp1 = dt*interp1(x,u,xp0,'pchip');
xp2 = dt*interp1(x,u,xp0+xp1/2,'pchip');
xp3 = dt*interp1(x,u,xp0+xp2/2,'pchip');
xp4 = dt*interp1(x,u,xp0+xp3,'pchip');

xp = xp0 + (xp1 + 2*xp2 + 2*xp3 + xp4)/6 + xoff;

% $$$ %*********************************************************************
% $$$ 
% $$$ function UI = interpolate(x, U, dx)
% $$$ 
% $$$ % Interpolate function U (velocity component u on grid) to
% $$$ % particle positions x. Uses Lagrangian interpolation of order
% $$$ % Iord (Iord = 1 ==> cubic interpolation).  See
% $$$ % Duran Ch. 6, for example.  ax below is x component
% $$$ % of what he calls alpha, the fractional grid position.
% $$$ 
% $$$ % nx = 2*(kmax+1)
% $$$ % dx = 2*pi/nx 
% $$$ 
% $$$ Iord = 1;
% $$$ bump = 10^(-10); % Prevent NaNs
% $$$ 
% $$$ nx = length(U);
% $$$ UI = zeros(size(x));  % size of output field = # of particles
% $$$ 
% $$$ for m=1:length(x)
% $$$         
% $$$     % get x on 0:2*pi-dx (remove periodic wrap-around), in units of
% $$$     % number of grid points.  e.g. notice that mod((2*pi-dx)/dx,nx)=nx-1
% $$$     xl = mod(x(m)/dx, nx);
% $$$         
% $$$     % get indeces of left/bottom grid point of cell containing
% $$$     % particle, min(i0) = 1
% $$$     i0 = 1 + floor(xl);
% $$$         
% $$$     % get fractional position within cell, 0 <= ax < 1
% $$$     ax = 1 + xl - i0;
% $$$     
% $$$     wx = ones(2*(Iord+1),1); 
% $$$     for i=-Iord:Iord+1
% $$$         for j=-Iord:Iord+1
% $$$             if (i~=j) 
% $$$                 wx(i+Iord+1) = wx(i+Iord+1)*(ax - j + bump)/(j - i);
% $$$             end 
% $$$         end 
% $$$     end 
% $$$     
% $$$     for i=-Iord:Iord+1
% $$$         for j=-Iord:Iord+1
% $$$             ig = 1 + mod( i0 + i - 1, nx );
% $$$             UI(m) = UI(m) + wx(i+Iord+1)*U(ig);
% $$$         end 
% $$$     end     
% $$$ end
% $$$       
% $$$ 
