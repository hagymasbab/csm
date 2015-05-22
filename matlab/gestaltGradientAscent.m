function gestaltGradientAscent(ge,data,batchSize,batchNum,stepNum,learningRate,priorSamples,randseed)
    if ge.B > 1 || ge.nullComponent
        error('not implemented');
    end

    setrandseed(randseed);
    N_all = size(data,1);
    data = reshape(data,N_all,ge.Dx);
    
    % initial condition
    % cc = eyes(ge.k,ge.Dv,0.001);    
    cc = randomCovariances(ge.k,ge.Dv);
    choles = cellchol(cc);
    
    loadSamples = false;
    state_sequence = {};
    state.estimated_components = cc;
    state_sequence{1} = state;
    
    batch_like = zeros(batchNum,stepNum+1);
    
    for batch = 1:batchNum
        printCounter(batch,'maxVal',batchNum,'stringVal','Batch');
        % sample a data batch
        img_indices = chooseKfromN(batchSize,N_all);
        X = data(img_indices,:);    
        batch_like(batch,1) = gestaltLogLikelihood2(ge,priorSamples,X,choles,'loadSamples',loadSamples);
        loadSamples = true;
        
        for step = 1:stepNum
            % calculate gradient
            grad = gestaltLogLikelihoodGradient(ge,priorSamples,X,choles,'loadSamples',loadSamples);            
            
%             gradmat = abs(cell2mat(grad));
%             max(gradmat(:))
            %pause
            % use the same set of samples at every iteration
            loadSamples = true;
            % update params 
            choles = celladd(choles,1,grad,learningRate);                                    
            
            viewImageSet(choles)
            % calculate likelihood on batch
            batch_like(batch,step+1) = gestaltLogLikelihood2(ge,priorSamples,X,choles,'loadSamples',loadSamples);        
            %pause
            save('gradasc_iter.mat','batch_like','ge','state_sequence','-v7.3');
        end
        % save stuff
        state.estimated_components = cholcell(choles);
        state_sequence{end+1} = state;
        save('gradasc_iter.mat','batch_like','ge','state_sequence','-v7.3');
    end
    
end