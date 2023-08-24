function [U, nablaU] = interpolate_U(background_flow1, background_flow2, alpha, x, h)
    xx = x(:, 1);
    yy = x(:, 2);
    
    U1(:,1) = interpolate(xx,yy,background_flow1.u,h,h);
    U1(:,2) = interpolate(xx,yy,background_flow1.v,h,h);
    U2(:,1) = interpolate(xx,yy,background_flow2.u,h,h);
    U2(:,2) = interpolate(xx,yy,background_flow2.v,h,h);
    
    nablaU1.u_x = interpolate(xx,yy,background_flow1.ux,h,h);
    nablaU1.u_y = interpolate(xx,yy,background_flow1.uy,h,h);
    nablaU1.v_x = interpolate(xx,yy,background_flow1.vx,h,h);
    nablaU1.v_y = interpolate(xx,yy,background_flow1.vy,h,h);
    nablaU2.u_x = interpolate(xx,yy,background_flow2.ux,h,h);
    nablaU2.u_y = interpolate(xx,yy,background_flow2.uy,h,h);
    nablaU2.v_x = interpolate(xx,yy,background_flow2.vx,h,h);
    nablaU2.v_y = interpolate(xx,yy,background_flow2.vy,h,h);
    
    U = (1 - alpha) * U1 + alpha * U2;
    nablaU.u_x = (1 - alpha) * nablaU1.u_x + alpha * nablaU2.u_x;
    nablaU.u_y = (1 - alpha) * nablaU1.u_y + alpha * nablaU2.u_y;
    nablaU.v_x = (1 - alpha) * nablaU1.v_x + alpha * nablaU2.v_x;
    nablaU.v_y = (1 - alpha) * nablaU1.v_y + alpha * nablaU2.v_y;
end