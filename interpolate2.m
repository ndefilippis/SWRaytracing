function FI = interpolate2(X_extend, Y_extend, F, Xq, Yq, L)
    num_extend_points = 4;    
    F_right = F(:,1:num_extend_points);
    F_left = F(:,end-num_extend_points+1:end);
    
    F_extend = [F_left, F, F_right];
    
    F_top = F_extend(end-num_extend_points+1:end,:);
    F_bottom = F_extend(1:num_extend_points,:);
    
    F_extend = [F_top; F_extend; F_bottom];
    
    %Xq_centered = mod(Xq + L/2, L) - L/2;
    %Yq_centered = mod(Yq + L/2, L) - L/2;
    Xq_centered = mod(Xq, L);
    Yq_centered = mod(Yq, L);
    
    FI = interp2(X_extend, Y_extend, F_extend', Xq_centered, Yq_centered, 'cubic');
end