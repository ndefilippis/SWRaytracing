

nx = 128;  
nt = 1000;
T = 50;
L = 2*pi;

dx = L/nx;
dt = T/nt;
x = linspace(0,L-dx,nx)-L/2;
t = linspace(0,T-dt,nt);

[x_,y_] = ndgrid(x,x);

np = 10;

psi = cos(x_).*cos(y_);
u = -cos(x_).*sin(y_);
v = sin(x_).*cos(y_);

x0 = (0:np-1)/np * L + 1e-7;
[foox,fooy] = ndgrid(x0,x0);  % positions of np^2 paticles at t=0
xp = foox(:);  
yp = fooy(:);  % use 1-d arrays
    
Xp(1).x = xp;
Xp(1).y = yp;
figure(11);
clf
set(gca,'fontsize',16)

clear M
for j=1:nt-1
    %[xp(:,j+1),yp(:,j+1)] = advect_particles(xp(:,j),yp(:,j),u,v,dt,dx); 
    [Xp(j+1)] = advect_particles(Xp(j),u,v,dt,dx); 
        
    clf;
    pcolor(x,x,psi'), shading interp, colorbar, axis image
    hold;
    plot(Xp(j+1).x,Xp(j+1).y,'k.','MarkerSize',10)
    drawnow
    M(j)=getframe;
end
xp=[Xp.x]; yp=[Xp.y];
