function r = symmetricDirichlet(alpha,dim,n)
    % take a sample from a dirichlet distribution
    a = alpha * ones(1,dim);
    r = gamrnd(repmat(a,n,1),1,n,dim);
    r = r ./ repmat(sum(r,2),1,dim);
end
