function plotGabor(f)
    dim = sqrt(size(f,2));
    im = reshape(f,dim,dim);
    [x,y] = meshgrid(1:1:dim,1:1:dim);
    pcolor(x,y,im);
    %shading('interp');
end