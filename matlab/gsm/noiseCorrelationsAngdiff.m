function noiseCorrelationsAngdiff(filters,priorcov,orients,sigma_x)
    z = 1;
    nFilt = length(orients);
    postCov = stableInverse(stableInverse(priorcov) + (z^2 / sigma_x^2) * (filters' * filters));
    %postCorr = corrcov(postCov);
    postCorr = corrcov(stableInverse(filters'*filters));
    %postCorr = corrcov(priorcov);
    
    corr_vs_orientdiff = [];
    
    for i = 1:nFilt
        for j = i+1:nFilt
            ordiff = abs(orients(i)-orients(j));
            corr_vs_orientdiff = [corr_vs_orientdiff; ordiff postCorr(i,j)];
        end
    end
    
    f = fit(corr_vs_orientdiff(:,1),corr_vs_orientdiff(:,2),'poly2');
    
    plot(f,corr_vs_orientdiff(:,1),corr_vs_orientdiff(:,2));
    ylim([-0.1 0.2])
        
end