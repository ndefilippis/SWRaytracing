function [U] = lsw1(Ui,f,C,time)

% Get initial spectral Uki from Ui
% For each k:
%   Create Vk = [U0k Upk Umk] and wk
%   Initial coefficients Ck = inv(Vk)*Uki
%   Uk(t) = C0*U0k + Cp*exp(i*wk*t)*Upk+ Cm*exp(-i*wk*t)*Umk
%
% U(t) = k2g(Uk)
% U(nx,nx,3)
% 
% H = [0    1i*f k*C; 
%     -1i*f 0    l*C; 
%      k*C  l*C  0  ];


nx = length(U,1);
if floor(log2(nx))~=log2(nx)
    error('make nx a power of 2');
end
kmax = (nx/2)-1;
kk = [0:kmax,-kmax-1:-1];
ll = [0:kmax,-kmax-1:-1];

%[k_,l_] = ndgrid(k,l);

Ui(:,3) = C*Ui(:,3);  % Trick to make operator Hermitian, a la Salmon 2.9
Uki = fft2(Ui);        % Get initial spectral fields

for ik=1:length(kk)
    for il=1:length(ll)
        [S(ik,il,:,:),w(ik,il)] = getVk(squeeze(Uki(ik,il,:)),f,C,kk(ik),ll(il));
    end
end

for t=time
    Uk = S(:,:,:,1) ...
       + S(:,:,:,2).*repmat(exp(-1i*w*t),[1 1 3]) ...
       + S(:,:,:,3).*repmat(exp( 1i*w*t),[1 1 3]);
    U(:,:,:,it) = ifft2(Uk);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,w] = getUk(Uki,f,C,k,l)

% Given column vector Uki = U(t=0,k,l) with Uki(:) = [u; v; C*h]
% get eigenfunctions and initial coefficients, to construct full
% solution at any time t:  
% Uk = Ck(1)*V0 + Ck(2)*Vp*exp(-1i*w*t) + Ck(3)*Vm*exp(1i*w*t);
% Define matrix S = [Ck(1)*V0, Ck(2)*Vp, Ck(3)*Vm]

K2 = k^2+l^2;
w = sqrt(f^2 + C^2*K2)     % Frequency array for each K


V0 = [-1i*l*C;      1i*k*C;      f];
Vp = [ w*k+1i*f*l;  w*l-1i*f*k;  C*K2];
Vm = [-w*k+1i*f*l; -w*l-1i*f*k;  C*K2];

V0 = V0/sqrt(V0'*V0);
Vp = Vp/sqrt(Vp'*Vp);
Vm = Vm/sqrt(Vm'*Vm);

V = [U0, Up, Um];  % checked that V' = inv(V)

Ck = V'*Uki;    % Initial coefficients

% solution for this wavenumber k,l and time t:
Uk = Ck(1)*V0 + Ck(2)*exp(-1i*w*t)*Vp + Ck(3)*exp(1i*w*t)*Vm;