function prodg = gprod(f,g)
    
    % Product of real parts and imaginary parts, loaded
    % respectively into real and imaginary part of output
        
    prodg = real(f).*real(g) + sqrt(-1)*imag(f).*imag(g);
