function [S,W,E] = getSk(Uki,f,C,k,l)

% [S,W,E] = getSk(Uki,f,C,k,l)
%
% Given column vector Uki = U(t=0,k,l) with Uki(:) = [u; v; C*h]
% get eigenfunctions, frequencies and initial coefficients to construct full
% solution at any time t:  
%   Uk = Ck(1)*V0 + Ck(2)*Vp*exp(-1i*W*t) + Ck(3)*Vm*exp(1i*W*t);
% Output matrix S = [Ck(1)*V0, Ck(2)*Vp, Ck(3)*Vm] 
% Also save energies E in each mode

K2 = k^2+l^2;
W = sqrt(f^2 + C^2*K2);     % Frequency array for each K

% Column vectors
V(:,1) = [-1i*l*C;      1i*k*C;      f];     % V0
V(:,2) = [ W*k+1i*f*l;  W*l-1i*f*k;  C*K2];  % V+
V(:,3) = [-W*k+1i*f*l; -W*l-1i*f*k;  C*K2];  % V-

% V is orthogonal, and so V' * V is diagonal, and entries are
% energies each mode.  Save these, then make V orthonormal
E = diag(V' * V);
for j = 1:3
    if E(j)~=0, V(:,j) = V(:,j)/E(j); end
end

Ck = V'*Uki;       % Initial coefficients

S = [Ck(1)*V(:,1), Ck(2)*V(:,2), Ck(3)*V(:,3)];

