function B = cc2B(cc)
    k = length(cc);
    Dv = size(cc{1},1);
    B = zeros(Dv,k);
    for kk = 1:k
        for i = 1:Dv
            for j = i+1:Dv
                B(i,k) = B(i,k) + cc{kk}(i,j);
                B(j,k) = B(j,k) + cc{kk}(i,j);
            end
        end
    end
end