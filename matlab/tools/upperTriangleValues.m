function B = upperTriangleValues(A)
%     A = A-diag(diag(A));
%     B = A(triu(true(size(A))));
%     B = B(B~=0);
    B = [];
    for i=1:size(A,1)
        for j=i+1:size(A,2)
            B = [B; A(i,j)];
        end
    end
end