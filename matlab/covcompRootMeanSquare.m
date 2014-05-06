function [mindiff,minperm] = covcompRootMeanSquare(cc1,cc2,minperm)
    % we should take the minimum over all possible permutations
    k = size(cc1,2);
    d = size(cc1{1},1);
    if isempty(minperm);
        permutations = perms(1:k);
    else
        permutations = minperm;
    end
    mindiff = Inf;
    for p=1:size(permutations,1)
        si = permutations(p,:);
        diff = 0;
        for j=1:k
            diff = diff + sum(sum((cc2{j}-cc1{si(1,j)})^2));           
        end
        diff = diff / (k*d^2);
        if diff < mindiff
            mindiff = diff;
            minperm = si;
        end
    end
end