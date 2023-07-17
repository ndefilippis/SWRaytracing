% addpath /home/qx344/SWcode/rsw
% addpath /home/qx344/SWcode/qgmat


nx = 256; % must be a power of 2

f = 3;   % = the deformation wavenumber k_d with Cg = 1
Cg = 1;

rng(14);
%----------------------------------------------------

frames = 100;
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx)' - L/2;
kmax = nx/2 - 1;  % 2^n - 1
nkx = 2*kmax+1;   % -kmax:kmax
dx = 2*pi/nx;     % L_domain = 2*pi
[x_,y_] = ndgrid(x,x);
tt = linspace(0,4*pi,frames);
%phi = rand*2*pi;
phi = randi([0, 1], [10,10])*2*pi;
Ufull = zeros(nx,nx,3,frames);
for T = 1:100
Uw = zeros(nx,nx,3);
for k = 1:10
    for l = 1:10
        K2 = k^2+l^2;
          % random phase
        theta = k*x_ + l*y_ + phi(k,l);
%         if (K2<= 5^2 && K2>3^2 && K2~=0)
%         if(k==4 && l==4)
%        if( K2==2*16 || K2==2*9 || K2==2*25)
            wp = sqrt(f^2+Cg^2*K2);
            theta = theta - wp*tt(T);
%             if rand>.5, wp = -wp; end  % half +, half -
            %aw = a0; %a0*f/wp;        
            etaw = cos(theta);
            uw = (k*wp*cos(theta)-l*f*sin(theta))/K2;
            vw = (l*wp*cos(theta)+k*f*sin(theta))/K2;
            Uw(:,:,1) = Uw(:,:,1) + uw;
            Uw(:,:,2) = Uw(:,:,2) + vw;
            Uw(:,:,3) = Uw(:,:,3) + etaw;
%         end
        % add a geostrophic component:
    end
end

Ufull(:,:,:,T) = Uw;
end


ag = 0.05;
aw = 0.1;

a = .25;
km = 5;

Psi = Cg^2/f*( sin(km*x_).*sin(km*y_)+a*cos(km*x_).*cos(km*y_) );
Psi = ag*Psi/max(max(abs(Psi(:,:,1))));

clear theta etaw uw vw 


Uinit = Ufull(:,:,:,100);
Uinit = aw*Uinit/max(max(abs(Uinit(:,:,1))));


%%
numsteps = 10000;
savestep = 25;
raXT = 30;
[U,t,ke,pe,Xp,hmov] = swkU_tc(Uinit,Psi,f,Cg,numsteps,savestep, x_, y_, km, a, ag, raXT);
save('U_0.mat','U');
save('t_0.mat','t');