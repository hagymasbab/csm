function comb = combinations(values,vecLength,actLength)
    
    if nargin < 3
        actLength = 0;
    end

    comb = [];
    if actLength < vecLength
        subcomb = combinations(values,vecLength,actLength+1);
        nSubComb = max(size(subcomb,1),1);
        for i = 1:length(values)
            newBlock = [values(i) * ones(nSubComb,1) subcomb];
            comb = [comb; newBlock];
        end
    end
end