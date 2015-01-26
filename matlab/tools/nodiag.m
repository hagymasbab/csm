function ndM = nodiag(M)    
    if ~iscell(M)
        M = {M};        
    end
    ndM = cell(size(M));
    for i = size(M,1)
        for j = size(M,2)
            ndM{i,j} = M{i,j} - diag(diag(M{i,j}));
        end
    end
    if length(ndM) == 1
        ndM = ndM{1};
    end
end