
nx = 128;  
nt = 1000;
T = 50;
L = 2*pi;

dx = L/nx;
dt = T/nt;
x = linspace(0,L-dx,nx)-L/2;
t = linspace(0,T-dt,nt);

u = (cos(x)).^2;

x0 = [0 2 4 6 14];
np = length(x0);
xp = zeros(np,nt);
xp(:,1) = x0;

for j=1:nt-1
    xp(:,j+1) = advect1d(xp(:,j),u,x,dt);
end

% dX/dt = u = cos^2(x)
% ==> int dX/u = tan(X)-tan(X0) = t 
% ==> X = atan(t+tan(X0)), with offsets of pi dep on X0
clear xpth
xoff=pi*floor((x0+pi/2)/pi);
for n=1:np
    xpth(n,:) = atan(t+tan(x0(n)))+xoff(n);
end

figure
plot(x,u)

figure
plot(t,xpth,'b')
hold
plot(t,xp,'r')

figure
plot(t,xpth-xp)




