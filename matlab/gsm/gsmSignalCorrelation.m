function sigCorr = gsmSignalCorrelation(nPatch,A,C,sigma_x)

    if ischar(A)
        load(A);
    end
    if ischar(C)
        load(C);
    end
    if ischar(sigma_x)
        load(sigma_x);
    end
        
    Dx = size(A,1);
    Dv = size(A,2);
    load(sprintf('patches_vanhateren_%d.mat',Dx));
    patches = patchDB(:,1:nPatch)';
    
    posteriorMeans = zeros(nPatch,Dv);
    
    for p=1:nPatch
        printCounter(p,'maxVal',nPatch,'stringVal','Patch');
        [mu_post,~,~,~] = gsmPosteriorV(patches(p,:)',A,C,sigma_x,2,2,30);
        posteriorMeans(p,:) = mu_post;
    end
    
    sigCorr = corr(posteriorMeans);
end