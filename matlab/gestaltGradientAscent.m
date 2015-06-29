function gestaltGradientAscent(ge,data,batchSize,batchNum,stepNum,varargin)
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'likeComp','none');
    addParameter(parser,'randseed','shuffle');      
    addParameter(parser,'priorSamples',100,@isnumeric); 
    addParameter(parser,'learningRate',0.1,@isnumeric); 
    addParameter(parser,'template',true,@islogical);
    addParameter(parser,'dampA',true,@islogical);
    addParameter(parser,'testLike',0);  
    addParameter(parser,'initCond','random');
    addParameter(parser,'sigmaSteps',0);
    addParameter(parser,'startWithSigma',false,@islogical);
    parse(parser,varargin{:});        
    params = parser.Results;  

    if ge.B > 1 || ge.nullComponent
        error('not implemented');
    end
    loadSamples = false;    
    setrandseed(params.randseed);
    N_all = size(data,1);
    data = reshape(data,N_all,ge.Dx);   
    
    verb = 0;
    if params.verbose >= 3
        verb = 1;
    end
    like_method = 'algebra';    
    
    t = true(ge.Dv);
    if params.template
        load('covariance_template_576.mat');
        %t = covarianceTemplate(ge.A,{'overlap','parallell'},{0.05,5});
    end   
    
    % initial condition    
    if strcmp(params.initCond,'random')
        cc = randomCovariances(ge.k,ge.Dv);
        choles = cellchol(cc);
        for i = 1:ge.k
            choles{i}(~t) = 0;
        end
        cc = cholcell(choles);
        state_sequence = {};
        state.estimated_components = cc;
        state.estimated_cholesky = choles;
        state_sequence{1} = state;
        full_like = [];
        batch_like = [];
        test_like = [];
        batch_indices = [];
        test_indices = [];
        
        if ischar(params.testLike) || params.testLike > 0
            if ischar(params.testLike)
                load(params.testLike)
            else
                test_indices = chooseKfromN(params.testLike,N_all);
            end
            X_test = data(test_indices,:);          
            ll = gestaltLogLikelihood2(ge,params.priorSamples,X_test,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);
            test_like = [test_like;ll];
            loadSamples = true;
        end

    else
        load(params.initCond)
        if isfield(state_sequence{end},'estimated_cholesky')
            choles = state_sequence{end}.estimated_cholesky;
            cc = cholcell(choles);
        else
            cc = state_sequence{end}.estimated_components;
            choles = cellchol(cc,'cheat',true);
        end
    end    
    ge.cc = cc;
    
    if params.dampA
        ge.A = ge.A + 0.05*eye(ge.Dv);
    end                                

    for batch = 1:batchNum
        if params.verbose == 1
            printCounter(batch,'maxVal',batchNum,'stringVal','Batch');
        elseif params.verbose >= 2
            fprintf('Batch %d/%d ',batchNum,batch);
        end
        % sample a data batch
        img_indices = chooseKfromN(batchSize,N_all);
        X = data(img_indices,:);        
        act_batch_like = [];
        if ~strcmp(params.likeComp,'none')
            ll = gestaltLogLikelihood2(ge,params.priorSamples,X,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);
            act_batch_like = [act_batch_like ll];
            loadSamples = true;
        end        
        
        for step = 1:stepNum
            if params.verbose == 2
                printCounter(step,'maxVal',stepNum,'stringVal','Step');
            elseif params.verbose >= 3
                fprintf('Step %d/%d ',stepNum,step);
            end
            % calculate gradient
            grad = gestaltLogLikelihoodGradient(ge,params.priorSamples,X,choles,'loadSamples',loadSamples,'method','scinot','template',t,'verbose',verb);
            
%             gradmat = abs(cell2mat(grad));            
%             max(gradmat(:))
%             min(gradmat(:))
            %pause
            % use the same set of samples at every iteration
            loadSamples = true;
            % update params 
            choles = celladd(choles,1,grad,params.learningRate);                                    
            
            % viewImageSet(choles)
            % calculate likelihood on batch
            if ~strcmp(params.likeComp,'none')
                ll = gestaltLogLikelihood2(ge,params.priorSamples,X,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);                 
                act_batch_like = [act_batch_like ll];
            end
            %viewImageSet(grad)
            %pause
            save('gradasc_iter.mat','batch_like','ge','state_sequence','-v7.3');
        end
        batch_like = [batch_like; act_batch_like];
        if strcmp(params.likeComp,'full')
            ll = gestaltLogLikelihood2(ge,params.priorSamples,data,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);            
            full_like = [full_like;ll];
        end
        if params.testLike > 0       
            ll = gestaltLogLikelihood2(ge,params.priorSamples,X_test,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);          
            test_like = [test_like;ll];
        end
        % save stuff
        state.estimated_components = cholcell(choles);
        state.estimated_cholesky = choles;
        state_sequence{end+1} = state;
        batch_indices = [batch_indices; img_indices];
        learningRate = params.learningRate;
        save('bin/gradasc_iter.mat','batch_like','full_like','test_like','ge','state_sequence','batch_indices','test_indices','learningRate','-v7.3');
    end    
end