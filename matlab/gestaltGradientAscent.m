function gestaltGradientAscent(ge,data,batchSize,batchNum,stepNum,varargin)
    parser = inputParser;   
    addParameter(parser,'verbose',0,@isnumeric);  
    addParameter(parser,'likeComp','none');
    addParameter(parser,'randseed','shuffle');      
    addParameter(parser,'priorSamples',100,@isnumeric); 
    addParameter(parser,'learningRate',0.1,@isnumeric); 
    addParameter(parser,'template',true,@islogical);
    addParameter(parser,'dampA',true,@islogical);
    addParameter(parser,'testLike',0,@isnumeric);  
    parse(parser,varargin{:});        
    params = parser.Results;  

    if ge.B > 1 || ge.nullComponent
        error('not implemented');
    end

    t = true(ge.Dv);
    if params.template
        load('covariance_template_576.mat');
        %t = covarianceTemplate(ge.A,{'overlap','parallell'},{0.05,5});
    end    
    
    if params.dampA
        ge.A = ge.A + 0.05*eye(ge.Dv);
    end
    
    setrandseed(params.randseed);
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
    test_like = zeros(batchNum+1,1);
    batch_indices = [];
    test_indices = [];
    verb = 0;
    if params.verbose >= 3
        verb = 1;
    end
    like_method = 'algebra';
    
    if params.testLike > 0
        test_indices = chooseKfromN(params.testLike,N_all);
        X_test = data(test_indices,:);          
        test_like(1) = gestaltLogLikelihood2(ge,params.priorSamples,X_test,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);
        loadSamples = true;
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
        if ~strcmp(params.likeComp,'none')
            ll = gestaltLogLikelihood2(ge,params.priorSamples,X,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);
            batch_like(batch,1) = ll;
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
            
            %gradmat = abs(cell2mat(grad));            
            %max(gradmat(:))
            %pause
            % use the same set of samples at every iteration
            loadSamples = true;
            % update params 
            choles = celladd(choles,1,grad,params.learningRate);                                    
            
            % viewImageSet(choles)
            % calculate likelihood on batch
            if ~strcmp(params.likeComp,'none')
                ll = gestaltLogLikelihood2(ge,params.priorSamples,X,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);
                batch_like(batch,step+1) = ll;        
            end
            %viewImageSet(grad)
            %pause
            save('gradasc_iter.mat','batch_like','ge','state_sequence','-v7.3');
        end
        if strcmp(params.likeComp,'full')
            ll = gestaltLogLikelihood2(ge,params.priorSamples,data,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);
            full_like(batch) = ll;
        end
        if params.testLike > 0       
            test_like(batch+1) = gestaltLogLikelihood2(ge,params.priorSamples,X_test,choles,'loadSamples',loadSamples,'verbose',verb,'method',like_method);            
        end
        % save stuff
        state.estimated_components = cholcell(choles);
        state_sequence{end+1} = state;
        batch_inidces = [batch_indices; img_indices];
        save('bin/gradasc_iter.mat','batch_like','full_like','test_like','ge','state_sequence','batch_indices','test_indices','-v7.3');
    end    
end