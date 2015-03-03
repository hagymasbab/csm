function [elements,x,y] = maxNElements(A,N,excluded)
    for e=1:size(excluded,1)
        A(excluded(e,1),excluded(e,2)) = -Inf;
    end
    [sortedValues,sortIndex] = sort(A(:),'descend');
    maxIndex = sortIndex(1:N);
    elements = A(maxIndex);
    [x,y] = ind2sub(size(A),maxIndex);
end