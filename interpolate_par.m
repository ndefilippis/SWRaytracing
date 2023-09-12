function FI = interpolate_par(x, y, F, dx, dy)

% Interpolate function F (velocity component u or v, on grid) to
% particle positions (x,y). Uses Lagrangian interpolation of order
% Iord (Iord = 1 ==> cubic interpolation).  See
% Duran Ch. 6, for example.  ax,ay below are x and y components
% of what he calls alpha, the fractional grid position.

% nx = 2*(kmax+1)
% dx = 2*pi/nx

Iord = 2;
bump = 10^(-10); % Prevent NaNs

[nx,ny] = size(F);
FI = zeros(size(x));  % size of output field = # of particles

function fi=compute_FI(x, y)
    xl = mod(x/dx, nx);
    yl = mod(y/dy, ny);
    
    % get indeces of left/bottom grid point of cell containing
    % particle, min(i0,j0) = 1,1
    i0 = 1 + floor(xl);
    j0 = 1 + floor(yl);
        
    % get fractional position within cell, 0 <= ax,ay < 1
    ax = 1 + xl - i0;
    ay = 1 + yl - j0;
    
    wx = ones(2*(Iord+1),1); wy = wx;
    for i=-Iord:Iord+1
        for j=-Iord:Iord+1
            if (i~=j) 
                wx(i+Iord+1) = wx(i+Iord+1)*(ax - j + bump)/(j - i);
                wy(i+Iord+1) = wy(i+Iord+1)*(ay - j + bump)/(j - i);
            end 
        end 
    end 
    
    fi = 0;
    for i=-Iord:Iord+1
        for j=-Iord:Iord+1
            ig = 1 + mod( i0 + i - 1, nx );
            jg = 1 + mod( j0 + j - 1, nx );
            fi = fi + wx(i+Iord+1)*wy(j+Iord+1)*F(ig,jg);
        end 
    end     
end

FI = arrayfun(@compute_FI, x, y);

end



      

