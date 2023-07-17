n = 6;
KMAX = 2^n-1;
NX = 2*(KMAX+1);
KMAXBIG = 3*(KMAX+1)/2-1;
NXBIG = 2*(KMAXBIG+1);

x = linspace(0,2*pi,NX+1); x=x(1:end-1);
xb = linspace(0,2*pi,NXBIG+1); xb=xb(1:end-1);
K = [0:KMAX -KMAX-1:-1];

fg = zeros(NX,3);
fg(:,1) = sin(2*x) -3*cos(3*x) +.1*sin(7*x);
fg(:,2) = 2*fg(:,1);
fg(:,3) = 3*fg(:,1);

fk = fft(fg)/NX;

uk = fk(1:KMAX+1,:);  % this is our test input

ukbig = zeros(KMAXBIG+1,3);
ukbig = [uk; zeros((KMAX+1)/2,3)];  % [0:KMAXBIG,3]
ukbigf = zeros(NXBIG,3);
ukbigf(1:KMAXBIG+1,:) = ukbig;   
ukbigf(KMAXBIG+3:end,:) = conj(ukbig(end:-1:2,:));
ugbig = real(ifft(ukbigf));
u2kbigf = fft(ugbig.^2);
u2k = u2kbigf(1:KMAX+1);

xsb = linspace(0,2*pi,10*NXBIG+1); xsb=xsb(1:end-1);
fsb = sin(2*xsb) -3*cos(3*xsb) +.1*sin(7*xsb);

figure(1)
clf
plot(x,fg(:,1),'o-')
hold
plot(xb,ugbig(:,1),'rs-')
plot(xsb,fsb,'k')
grid

