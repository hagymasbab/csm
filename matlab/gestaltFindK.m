function gestaltFindK(Dx,filterset,krange,dataset,trainingSize,batchSize,testSize,emSteps,samples)
    
    dataSize = trainingSize + testSize;
    nk = length(krange);

    if isnumeric(dataset)
        ge = gestaltCreate('temp','Dx',Dx,'k',dataset,'B',1,'N',dataSize,'filters',filterset, ...
            'obsVar',0.1,'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'nullComponent',false,'generateComponents',true,'generateData',true);
        X = reshape(ge.X,dataSize,Dx);
        clear ge;
        sciLike = false;
    else
        datafile = sprintf('patches_%s_%d.mat',dataset,Dx);
        load(datafile);
        X = patchDB(:,1:dataSize)';
        clear patchDB;
        sciLike = true;
    end
    
    X_train = X(1:trainingSize,:);
    X_test = X(trainingSize+1:end,:);
    
    test_likelihoods = zeros(nk,1);
    
    for i = 1:nk
        act_k = krange(i);
        ge = gestaltCreate('temp','Dx',Dx,'k',act_k,'B',1,'N',batchSize,'filters',filterset, ...
            'obsVar',1,'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'nullComponent',false,'generateComponents',false,'generateData',false);
        
        fprintf('Actual k: %d ',act_k);
        gestaltEM(ge,X_train,batchSize,emSteps,samples,'shuffle','syntheticData',false,'burnin',20,'savingCode',i);
        load(sprintf('cc_%d.mat',i));
        ge.cc = cc_saved;
        test_likelihoods(i) = gestaltLogLikelihood(ge,100,X_test,'scientific',sciLike);
        save('findk_likes.mat','test_likelihoods','krange','ge');
        fprintf('\b test set log-likelihood: %.4f\n',test_likelihoods(i));
    end
end