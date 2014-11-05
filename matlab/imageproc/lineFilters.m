function A = lineFilters(imSize)
    Dx = imSize^2;
    A = zeros(Dx,Dx);
    template = [0.5 0.9 0.5];    
    length = 3;
    contrast = 2;
    lowval = 1 / (2*length+contrast);
    highval = contrast * lowval;
    for i=1:imSize
        for j=1:imSize
            minc = max(1,j-length);
            maxc = min(imSize,j+length);
            rowidx = (i-1)*imSize + j;
            A(rowidx,(i-1)*imSize + minc : (i-1)*imSize + maxc) = lowval;
            A(rowidx,rowidx) = highval;
        end
    end    
end