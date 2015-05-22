function cellsum = celladd(c1,s1,c2,s2)
    k = size(c1,2);
    cellsum = cell(1,k);
    for i=1:k
        cellsum{i} = c1{i} * s1 + c2{i} * s2;
    end
end