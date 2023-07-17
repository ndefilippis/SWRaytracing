function [Pout] = step_packet_xka(P,U,GradU,H,C0,f,dx,dy,dt)
    
% [Pout] = step_packet(P,U,GradU,H,C0,f,dx,dy,dt)
%
% Advects one wavepacket for time dt with RK4 scheme, solving
%
% dX/dt = u + omega_k 
% dY/dt = v + omega_l
% dk/dt = -u_x k - v_x l - omega_x
% dl/dt = -u_y k - v_y l - omega_y
% da/dt = -a * div(C)
% 
% (neglecting dependence of U, omega on t)  
%
% The 2D velocity field and its gradients are
% specified on grid with spacing dx, dy, and passed as structures:
%
% u = U.u(:,:)
% v = U.v(:,:) 
% u_x = GradU.u_x(:,:)
% u_y = GradU.u_y(:,:)
% v_x = GradU.v_x(:,:)
% v_y = GradU.v_y(:,:)
%
% Packet position x, y, and wavenumbers k, l are specified by
% structure P with fields P.x, P.y, P.k, P.l   
%
% Note that for RSW with Doppler and refraction, any strain in the
% flow makes all three frequency roots complex, but vorticity alone
% is okay.  Then
%
% omega^2 = f_eff^2 + K^2*c0^2, f_eff = f + zeta/2
%
% ... but this is for later

% Get group velocity for this packet's k, l

[C,~,~,divC,gradomega] = cg_sw(P.k,P.l,C0,f,U,H); 

% Step position of packet forward to new x, y

x1 = dt*interpolate(P.x,P.y,U.u+C.x,dx,dy);
y1 = dt*interpolate(P.x,P.y,U.v+C.y,dx,dy);

x2 = dt*interpolate(P.x+x1/2,P.y+y1/2,U.u+C.x,dx,dy);
y2 = dt*interpolate(P.x+x1/2,P.y+y1/2,U.v+C.y,dx,dy);

x3 = dt*interpolate(P.x+x2/2,P.y+y2/2,U.u+C.x,dx,dy);
y3 = dt*interpolate(P.x+x2/2,P.y+y2/2,U.v+C.y,dx,dy);

x4 = dt*interpolate(P.x+x3,P.y+y3,U.u+C.x,dx,dy);
y4 = dt*interpolate(P.x+x3,P.y+y3,U.v+C.y,dx,dy);

Pout.x = P.x + (x1 + 2*x2 + 2*x3 + x4)/6;
Pout.y = P.y + (y1 + 2*y2 + 2*y3 + y4)/6; 

% Now interpolate gradients to new positions

u_xi = interpolate(Pout.x,Pout.y,GradU.u_x,dx,dy);
u_yi = interpolate(Pout.x,Pout.y,GradU.u_y,dx,dy);
v_xi = interpolate(Pout.x,Pout.y,GradU.v_x,dx,dy);
v_yi = interpolate(Pout.x,Pout.y,GradU.v_y,dx,dy);
omega_xi = interpolate(Pout.x,Pout.y,gradomega.x,dx,dy);
omega_yi = interpolate(Pout.x,Pout.y,gradomega.y,dx,dy);
divCi = interpolate(Pout.x,Pout.y,divC,dx,dy);

% Step wavenumbers forward

k1 = dt*(-u_xi*P.k - v_xi*P.l - omega_xi);
l1 = dt*(-u_yi*P.k - v_yi*P.l - omega_yi);

k2 = dt*(-u_xi*(P.k+k1/2) - v_xi*(P.l+l1/2)- omega_xi);
l2 = dt*(-u_yi*(P.k+k1/2) - v_yi*(P.l+l1/2)- omega_yi);

k3 = dt*(-u_xi*(P.k+k2/2) - v_xi*(P.l+l2/2)- omega_xi);
l3 = dt*(-u_yi*(P.k+k2/2) - v_yi*(P.l+l2/2)- omega_yi);

k4 = dt*(-u_xi*(P.k+k3) - v_xi*(P.l+l3)- omega_xi);
l4 = dt*(-u_yi*(P.k+k3) - v_yi*(P.l+l3)- omega_yi);

Pout.k = P.k + (k1 + 2*k2 + 2*k3 + k4)/6;
Pout.l = P.l + (l1 + 2*l2 + 2*l3 + l4)/6;

% Step action

a1 = dt*(-P.a*divCi);
a2 = dt*(-(P.a+a1/2)*divCi);
a3 = dt*(-(P.a+a2/2)*divCi);
a4 = dt*(-(P.a+a3)*divCi);

Pout.a = P.a + (a1 + 2*a2 + 2*a3 + a4)/6;


