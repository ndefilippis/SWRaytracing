function y=apply_3d(x, f)
    first = f(x(:,:,1));
    y = zeros([size(first), size(x, 3)]);
    for i=1:size(x, 3)
        y(:,:,i) = f(x(:,:,i));
    end
end