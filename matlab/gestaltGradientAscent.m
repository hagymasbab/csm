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
    
    for batch = 1:batchNum
        % sample a data batch
        img_indices = chooseKfromN(batchSize,N_all);
        X = data(img_indices,:);    
        ll = gestaltLogLikelihood2(ge,priorSamples,X,choles,'loadSamples',loadSamples)
        loadSamples = true;
        
        for step = 1:stepNum
            % calculate gradient
            grad = gestaltLogLikelihoodGradient(ge,priorSamples,X,choles,'loadSamples',loadSamples);            
            
            gradmat = abs(cell2mat(grad));
            max(gradmat(:))
            %pause
            % use the same set of samples at every iteration
            loadSamples = true;
            % update params 
            choles = celladd(choles,1,grad,learningRate);
            % calculate likelihood on batch
            ll = gestaltLogLikelihood2(ge,priorSamples,X,choles,'loadSamples',loadSamples)
            viewImageSet(choles)
            pause
        end
    end
    
end