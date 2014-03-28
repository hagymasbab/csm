function cc = eyes(k,dim)
    cc = cell(1,k);
    for i=1:k
        cc{i} = eye(dim);
    end
end