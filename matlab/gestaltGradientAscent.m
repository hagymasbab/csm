function gestaltGradientAscent(ge,template,data,batchSize,batchNum,stepNum,learningRate,priorSamples,randseed,like_comp,verbose)
    if ge.B > 1 || ge.nullComponent
        error('not implemented');
    end

    t = true(ge.Dv);
    if template
        t = covarianceTemplate(ge.A,{'overlap','parallell'},{0.05,5});
    end
    
    
    setrandseed(randseed);
    N_all = size(data,1);
    data = reshape(data,N_all,ge.Dx);
    
    % initial condition
    % cc = eyes(ge.k,ge.Dv,0.001);    
    cc = randomCovariances(ge.k,ge.Dv);
    choles = cellchol(cc);
    for i = 1:ge.k
        choles{i}(~t) = 0;
    end
    cc = cholcell(choles);
    ge.cc = cc;
    
    loadSamples = false;
    state_sequence = {};
    state.estimated_components = cc;
    state_sequence{1} = state;
    
    batch_like = zeros(batchNum,stepNum+1);
    full_like = zeros(batchNum,1);
    

    for batch = 1:batchNum
        if verbose == 1
            printCounter(batch,'maxVal',batchNum,'stringVal','Batch');
        elseif verbose == 2
            fprintf('Batch %d/%d ',batchNum,batch);
        end
        % sample a data batch
        img_indices = chooseKfromN(batchSize,N_all);
        X = data(img_indices,:);           
        if ~strcmp(like_comp,'none')
            ll = gestaltLogLikelihood2(ge,priorSamples,X,choles,'loadSamples',loadSamples);            
            batch_like(batch,1) = ll;
            loadSamples = true;
        end        
        
        for step = 1:stepNum
            if verbose == 2
                printCounter(step,'maxVal',stepNum,'stringVal','Step');
            end
            % calculate gradient
            grad = gestaltLogLikelihoodGradient(ge,priorSamples,X,choles,'loadSamples',loadSamples,'method','scinot','template',t);
            
            %gradmat = abs(cell2mat(grad));            
            %max(gradmat(:))
            %pause
            % use the same set of samples at every iteration
            loadSamples = true;
            % update params 
            choles = celladd(choles,1,grad,learningRate);                                    
            
            % viewImageSet(choles)
            % calculate likelihood on batch
            if ~strcmp(like_comp,'none')
                ll = gestaltLogLikelihood2(ge,priorSamples,X,choles,'loadSamples',loadSamples);
                batch_like(batch,step+1) = ll;        
            end
            %viewImageSet(grad)
            %pause
            save('gradasc_iter.mat','batch_like','ge','state_sequence','-v7.3');
        end
        if strcmp(like_comp,'full')
            ll = gestaltLogLikelihood2(ge,priorSamples,data,choles,'loadSamples',loadSamples);
            full_like(batch) = ll;
        end
        % save stuff
        state.estimated_components = cholcell(choles);
        state_sequence{end+1} = state;
        save('gradasc_iter.mat','batch_like','ge','state_sequence','-v7.3');
    end
    
end