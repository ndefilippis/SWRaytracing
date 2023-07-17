function [h,u,v] = getwavemodes(Kv,A,sgn,x,y,gam,ep)

%  [h,u,v] = getwavemodes(Kv,A,sgn,x,y,gam,ep)
%  Get real-space RSW u,v,h for linear superposition of Nw waves, with 
%
%  h = Re{ Sum_j A(j) * exp(i*(kx(j)*x+ky(j)*y))}  
%
%  and dispersion relation 
%
%  omega(j) = sgn(j)*sqrt(1+gam*ep*K^2),  sgn = +1 or -1.
% 
%  Inputs: Kv(:,1) = kx(:), Kv(:,2) = ky(:).  
%  Number Nw of waves = size(Kv,1). 

h = zeros(size(x));
u = h; v = h;

for j=1:size(Kv,1)
    
    kx = Kv(j,1); 
    ky = Kv(j,2);
    K2 = kx^2+ky^2; 
    omega = sgn(j)*sqrt(1+gam*ep*K2);

    ht = A(j)*exp(1i*(kx*x+ky*y));
    h = h + ht;
    u = u + (omega*kx + 1i*ky)/(ep*gam*K2)*ht;
    v = v + (omega*ky - 1i*kx)/(ep*gam*K2)*ht;

end
    
h = (h + conj(h))/2; 
u = (u + conj(u))/2; 
v = (v + conj(v))/2;
 
