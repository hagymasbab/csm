function cholesky = cellchol(cc)
    cholesky = cc;
    for j=1:size(cc,2)
        cholesky{j} = chol(cholesky{j});
    end    
end