function ge = gestaltTestLearning(Dx,k,N,emBatchSize,nTrials,nSteps,nSamples,dataset,karklin,like_full)    
            
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
    else
        datafile = sprintf('patches_%s_%d.mat',dataset,Dx);
        load(datafile);
        X = reshape(patchDB(:,1:N)',N,1,Dx);
        synDat = false;
    end    
        
    ll = zeros(nTrials,nSteps+1);
    comps = cell(nTrials,nSteps+1,ge.k);
    
    if like_full
        like_comp = 'full';
    else
        like_comp = 'batch';
    end
    
    for t=1:nTrials
        printCounter(t,'stringVal','Trial','maxVal',nTrials);        
        gestaltEM(ge,X,emBatchSize,nSteps,nSamples,'shuffle','syntheticData',synDat,'burnin',20,'verbose',0,'cctComponents',karklin,'computeLikelihood',like_comp);
        load iter;
        for s=1:nSteps+1
            if like_full
                ll(t,s) = state_sequence{s}.full_loglike;
            else
                if s==1
                    ll(t,s) = state_sequence{s+1}.batch_loglike_pre;
                else
                    ll(t,s) = state_sequence{s}.batch_loglike;
                end
            end
            for kk=1:ge.k
                comps{t,s,kk} = state_sequence{s}.estimated_components{kk};
            end
        end
        
        save('comps.mat','comps');
        save('loglikes.mat','ll','like_comp');
    end
    
end