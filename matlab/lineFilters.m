function A = lineFilters(imSize)
    Dx = imSize^2;
    A = zeros(Dx,Dx);
    template = [0.5 0.9 0.5];
    length = floor(template/2);
    for i=1:imSize
        for j=1:imSize
            minc = max(1,j-length);
            maxc = min(imSize,j+length);
            rowidx = (i-1)*imSize + j;
            A(rowidx,(i-1)*imSize + minc : (i-1)*imSize + maxc) = 0.5;
            A(rowidx,rowidx) = 0.9;
        end
    end
end