nx = 256; % must be a power of 2
L = 2*pi;
x = linspace(0,L*(nx-1)/nx,nx)' - L/2;
kmax = nx/2 - 1;  % 2^n - 1
nkx = 2*kmax+1;   % -kmax:kmax
dx = 2*pi/nx;     % L_domain = 2*pi

[x_,y_] = ndgrid(x,x);

f = 3;   % = the deformation wavenumber k_d with Cg = 1
Cg = 1;

ag = 0.3;
aw = 0.1;
kmi = 5;
kl = 10;
ku = 13;

% Now initialize separately a spectrum of waves with random phases, and
% random frequency branches, and a narrow-band geostrophic flow with 
% random phases.  

Uw = zeros(nx,nx,3);
Ug = Uw;
for k = -kmax:kmax
    for l = 0:kmax
        K2 = k^2+l^2;
        phi = rand*2*pi;  % random phase
        theta = k*x_ + l*y_ + phi;
        if (K2<= kmi^2 & K2~=0)
            wp = sqrt(f^2+Cg^2*K2);
            if rand>.5, wp = -wp; end  % half +, half -
            %aw = a0; %a0*f/wp;        
            etaw = cos(theta);
            uw = (k*wp*cos(theta)-l*f*sin(theta))/K2;
            vw = (l*wp*cos(theta)+k*f*sin(theta))/K2;
            Uw(:,:,1) = Uw(:,:,1) + uw;
            Uw(:,:,2) = Uw(:,:,2) + vw;
            Uw(:,:,3) = Uw(:,:,3) + etaw;
        end
        % add a geostrophic component:
        if (K2>kl^2 & K2<=ku^2)
            etag = cos(theta);
            ug =  Cg^2/f*l*sin(theta);
            vg = -Cg^2/f*k*sin(theta);
            Ug(:,:,1) = Ug(:,:,1) + ug;
            Ug(:,:,2) = Ug(:,:,2) + vg;
            Ug(:,:,3) = Ug(:,:,3) + etag;
        end
    end
end
Ug = ag*Ug/max(max(abs(Ug(:,:,1))));
Uw = aw*Uw/max(max(abs(Uw(:,:,1))));
Uin = Uw + Ug;

clear Ug Uw theta etag etaw uw vw ug vg x_ y_

np = 20;
Xpin.x = .2*randn(np,1);
Xpin.y = .2*randn(np,1);

[U,t,ke,pe,Xp,hmov] = swk(Uin,f,Cg,5000,100,Xpin);

Xa = mean(sqrt(Xp.x.^2+Xp.y.^2));
plot(t,Xa)

off=50;
q = getswpv(U(:,:,:,1:off:end),f);
delta = getdiv(U(:,:,:,1:off:end));

% KE

for j=1:size(U,4)
    uk = g2k(squeeze(U(:,:,1,j)));  
    vk = g2k(squeeze(U(:,:,2,j)));
    kes(:,j) = iso_spectra(uk.*conj(uk)+vk.*conj(vk));
end

% Get wave and geostrophic energies and spectra
[KEw,PEw,KEg,PEg,kes,pes,kesw,pesw,kesg,pesg,Um,etawk] = wavevortdecomp(U(:,:,:,end),f,Cg);

figure
loglog(kesw)
hold
loglog(kesg)
grid
set(gca,'fontsize',14)
xlabel('Wavenumber')
legend('KE_w','KE_w')
title('KE spectra')





% To do:
% Compute full Stokes drift from stable spectrum of waves:
% us = sum_k a_k^2*w_k/(2*k);
% where a_k are coefficients for FFT of eta (maybe with factor of 2)
% Full k-omega spectrum of u...  is it all still exact waves?
% time-avg of low-freq parts of stokes?  
% Work this all out for spectrum of waves in 2D SW.
% See if exact solution in 1D can be gotten in Fourier space, or at
% least understood that way, since single wave ends up as stable
% spectrum with a few wave modes. 
% Show shock as higher-order asymptotic solution...

