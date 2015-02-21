function [elements,x,y] = maxNElements(A,N)
    [sortedValues,sortIndex] = sort(A(:),'descend');
    maxIndex = sortIndex(1:N);
    elements = A(maxIndex);
    [x,y] = ind2sub(size(A),maxIndex);
end