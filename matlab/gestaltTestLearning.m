function gestaltTestLearning(Dx,k,N,emBatchSize,nTrials,nSteps,nSamples,filterset,dataset)    
            
    filterfile = sprintf('filters_%s_%d.mat',filterset,Dx);
       
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',N,'filters','filters_OF_64.mat', ...
        'obsVar',0.1,'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'nullComponent',false,'generateComponents',true,'generateData',true);

    if strcmp(dataset,'synthetic')
        X = ge.X;
        synDat = true;
    else
        datafile = sprintf('patches_%s_%d.mat',dataset,Dx);
        load(datafile);
        X = patchDB(:,1:N);
        synDat = false;
    end    
        
    ll = zeros(nTrials,nSteps+1);
    
    for t=1:nTrials
        printCounter(t,'stringVal','Trial','maxVal',nTrials);
        gestaltEM(ge,X,emBatchSize,nSteps,nSamples,'shuffle','syntheticData',synDat,'burnin',20,'verbose',0);
        load iter;
        for s=1:nSteps+1
            ll(t,s) = state_sequence{s}.loglike;
        end
        save('loglikes.mat','ll');
    end
    
end