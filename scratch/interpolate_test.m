L = 2*pi;
nx = 7;
x = linspace(-L/2, L/2, nx);
dx = L / (nx - 1);

x_interp = linspace(-L/2 - 4*dx, L/2 + 4*dx, nx + 8);

[X, Y] = meshgrid(x, x);
[X_interp, Y_interp] = meshgrid(x_interp, x_interp);

x = linspace(-L, L, 256);
y = linspace(-L/2, L/2, 256);

[test_X, test_Y] = meshgrid(x, y);
Z = interpolate2(X_interp, Y_interp, cos(X).*sin(Y), test_X, test_Y, L);
subplot(2, 1, 1);
contourf(X, Y, cos(X).*sin(Y));
colorbar();
subplot(2, 1, 2);
contourf(test_X, test_Y, reshape(Z, size(test_X)));
colorbar();
      

