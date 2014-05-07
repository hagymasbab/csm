function likelihood = gestaltLikelihood(ge,L)
    % approximated, up to a constant
    % get L samples from a k-dimensional symmetric dirichelet prior for g
    G = symmetricDirichlet(ge.sparsity,ge.k,L);
    likelihood = 0;
    for s=1:L        
        Cv = componentSum(G(s,:)',ge.cc);
        C = ge.obsVar * eye(ge.Dx) + ge.A * Cv * ge.A';
        dataProb = 0;
        for n=1:ge.N
            batchProb = 0;
            for b=1:ge.B                
                x = squeeze(ge.X(n,b,:));                                
                batchProb = batchProb + log(mvnpdf(x,zeros(size(x)),C));
                %pause
            end
            batchProb = exp(batchProb);
            dataProb = dataProb + batchProb;
        end        
        likelihood = likelihood + exp(dataProb);
    end
    likelihood = likelihood / L;
end