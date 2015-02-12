function gestaltTestLearning(Dx,k,N,emBatchSize,nTrials,nSteps,nSamples,dataset,karklin)    
            
    %filterfile = sprintf('filters_%s_%d.mat',filterset,Dx);
    
    if strcmp(dataset,'synthetic')
        genComps = true;
        genData = true;
        filterset = 'gabor_4or';
    else
        genComps = false;
        genData = false;
        filterset = 'OF';
    end
    
    ge = gestaltCreate('temp','Dx',Dx,'k',k,'B',1,'N',N,'filters',filterset, ...
        'obsVar',0.1,'g_shape',2,'g_scale',2,'z_shape',2,'z_scale',2,'nullComponent',false, ...
        'generateComponents',genComps,'generateData',genData,'componentShape','oriented-gabors');

    if strcmp(dataset,'synthetic')
        X = ge.X;
        synDat = true;
        comps = cell(nTrials,nSteps+1,ge.k);
    else
        datafile = sprintf('patches_%s_%d.mat',dataset,Dx);
        load(datafile);
        X = reshape(patchDB(:,1:N)',N,1,Dx);
        synDat = false;
    end    
        
    ll = zeros(nTrials,nSteps+1);
    
    for t=1:nTrials
        printCounter(t,'stringVal','Trial','maxVal',nTrials);
        gestaltEM(ge,X,emBatchSize,nSteps,nSamples,'shuffle','syntheticData',synDat,'burnin',20,'verbose',0,'cctComponents',karklin);
        load iter;
        for s=1:nSteps+1
            ll(t,s) = state_sequence{s}.loglike;
            if strcmp(dataset,'synthetic')
                %state_sequence{s}.estimated_components
                for kk=1:ge.k
                    comps{t,s,kk} = state_sequence{s}.estimated_components{kk};
                end
            end
        end
        
        if strcmp(dataset,'synthetic')
            save('comps.mat','comps');
        end
        save('loglikes.mat','ll');
    end
    
end