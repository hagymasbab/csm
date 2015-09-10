function B = upperTriangleValues(A,inclDiag)
%     A = A-diag(diag(A));
%     B = A(triu(true(size(A))));
%     B = B(B~=0);
    if nargin < 2
        inclDiag = false;
    end
    B = [];
    for i=1:size(A,1)
        if inclDiag
            begin = i;
        else
            begin = i+1;
        end
        for j=begin:size(A,2)
            B = [B; A(i,j)];
        end
    end
end