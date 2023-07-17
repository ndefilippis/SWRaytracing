function [U] = lsw(Ui,f,C,time)

% [U] = lsw(Ui,f,C,time)
% 
% Get Linear RSW solution at times in vector time, as follows
% 
% Multiply h = Ui(:,:,3) by C, so that propagation matrix is
% Hermitian:
% L = [0    1i*f k*C; 
%     -1i*f 0    l*C; 
%      k*C  l*C  0  ];
% with i * dU/dt = L * U, or L * U = w * U
% Get initial spectral Uki from Ui
% For each k:
%   Get eigenvalues w0 = 0, wp = W, wm = -W
%     where W = sqrt(f^2 + C^2*K^2),
%     and eigenvectors V = [V0 Vp Vm], where b/c L is Hermitian,
%     inv(V) = V'.
%   Initial coefficients Ck = V'*Uki
%   Full solution is
%   Uk(t) = Ck(1)*V0 + Ck(2)*Vp*exp(i*W*t)+ Ck(3)*Vm*exp(-i*W*t)
%   So create S = [Ck(1)*V0, Ck(2)*Vp, Ck(3)*Vm] , then get Uk(t)
%   for any t, and collect grid space output for all requested t = time
% U(t) = ifft2(Uk) which has size U(nx,nx,3,nt)
% 

nt = length(time);
nx = size(Ui,1);
if floor(log2(nx))~=log2(nx)
    error('make nx a power of 2');
end
kmax = (nx/2)-1;
kk = [0:kmax,-kmax-1:-1];
ll = [0:kmax,-kmax-1:-1];

Ui(:,3) = C*Ui(:,3);  % Trick to make operator Hermitian, a la Salmon 2.9
Uki = fft2(Ui);       % Get initial spectral fields

S = zeros(nx,nx,3,3);
W = zeros(nx,nx);
for ik=1:length(kk)
    for il=1:length(ll)
        [S(ik,il,:,:),W(ik,il)] = getSk(squeeze(Uki(ik,il,:)),f,C,kk(ik),ll(il));
    end
end

U = zeros(nx,nx,3,nt);
for it=1:nt
    t = time(it);
    Uk = S(:,:,:,1) ...
       + S(:,:,:,2).*repmat(exp(-1i*W*t),[1 1 3]) ...
       + S(:,:,:,3).*repmat(exp( 1i*W*t),[1 1 3]);
    U(:,:,:,it) = real(ifft2(Uk));
    U(:,:,3,it) = U(:,:,3,it)/C;
end
