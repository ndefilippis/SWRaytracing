n = 6;
KMAX = 2^n-1;
NX = 2*(KMAX+1);

x = linspace(0,2*pi,NX+1); 
x = x(1:end-1);

etahat = .05;
Bu = 1;
Ro = 1;
k = 4;
w = sqrt(1 + Bu*k^2);
c = w/k;
nu = 1e-16;

Uin = zeros(NX,3);
Uin(:,1) = c*etahat * cos(k*x);
Uin(:,2) = etahat/k * sin(k*x);
Uin(:,3) = etahat * cos(k*x);

[U,KE,PE,t] = sw1rk3nu(Uin,Ro,Bu,nu,1000,10);

figure
clear M
for j=1:size(U,3)
    plot(x,U(:,3,j))
    M(j) = getframe;
end

