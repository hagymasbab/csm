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
    
    full_like = [];
    batch_like = [];
    test_like = [];
    batch_indices = [];
    test_indices = [];
    verb = 0;
    if params.verbose >= 3
        verb = 1;
    end
    like_method = 'algebra';
    
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
        state_sequence{end+1} = state;
        batch_indices = [batch_indices; img_indices];
        save('bin/gradasc_iter.mat','batch_like','full_like','test_like','ge','state_sequence','batch_indices','test_indices','-v7.3');
    end    
end

%
% NEW STUFF
% 
% function gestaltGradientAscent(ge,data,batchSize,stepNum,varargin)
%     parser = inputParser;   
%     addParameter(parser,'verbose',0,@isnumeric);  
%     addParameter(parser,'likeComp','none');
%     addParameter(parser,'randseed','shuffle');      
%     addParameter(parser,'priorSamples',100,@isnumeric); 
%     addParameter(parser,'learningRate',0,@isnumeric); 
%     addParameter(parser,'template',true,@islogical);
%     addParameter(parser,'dampA',true,@islogical);
%     addParameter(parser,'testLike',0);  
%     addParameter(parser,'initCond','random');
%     addParameter(parser,'sigmaSteps',0);
%     addParameter(parser,'initSigma',0.5,@isnumeric); % only used if we start a new run
%     addParameter(parser,'startWithSigma',false,@islogical);
%     addParameter(parser,'synthetic',false,@islogical);
%     addParameter(parser,'likeMethod','intuition');
%     parse(parser,varargin{:});        
%     params = parser.Results;  
%     
%     if ge.B > 1 || ge.nullComponent
%         error('not implemented');
%     end
%     loadSamples = false;    
%     setrandseed(params.randseed);
%     N_all = size(data,1);
%     data = reshape(data,N_all,ge.Dx);   
%     batchNum = 1000;
%     if params.dampA
%         ge.A = ge.A + 0.05*eye(ge.Dv);
%     end 
%     
%     trueCC = {};
%     trueSigma = 0;
%     if params.synthetic
%         trueCC = ge.cc;
%         trueSigma = ge.obsVar;
%     end
%     
%     verb = 0;
%     if params.verbose >= 3
%         verb = 1;
%     end
%     likefunc = @(X_param,choles_param,sigma_param) gestaltLogLikelihood2(ge,params.priorSamples,X_param,choles_param,'loadSamples',loadSamples,'verbose',verb,'method',params.likeMethod,'sigma',sigma_param);
%     %likefunc = @(X_param,choles_param,sigma_param) gestaltLogLikelihood2(ge,params.priorSamples,X_param,choles_param,'loadSamples',loadSamples,'verbose',verb,'method',params.likeMethod);
%     
%     t = true(ge.Dv);
%     if params.template
%         load(sprintf('covariance_template_%d.mat',ge.Dv));
%         %t = covarianceTemplate(ge.A,{'overlap','parallell'},{0.05,5});
%     end   
%     
%     if params.learningRate == 0
%         % we don't take into account L and Dv
%         params.learningRate = 5 / batchSize;
%     end        
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%% SET INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     if strcmp(params.initCond,'random')
%         
%         ge.obsVar = params.initSigma;
%         
%         cc = randomCovariances(ge.k,ge.Dv);
%         choles = cellchol(cc);
%         for i = 1:ge.k
%             choles{i}(~t) = 0;
%         end
%         cc = cholcell(choles);
%         ge.cc = cc;   
%         state_sequence = {};
%         state.estimated_components = cc;
%         state.estimated_cholesky = choles;
%         state.sigma_x = ge.obsVar;
%         state.learningRate = params.learningRate;
%         state_sequence{1} = state;
%         full_like = [];
%         batch_like = [];
%         test_like = [];
%         batch_indices = [];
%         test_indices = [];
%         sigset_indices = [];
%         
%         if ischar(params.testLike) || params.testLike > 0
%             if ischar(params.testLike)
%                 load(params.testLike)
%             else
%                 test_indices = chooseKfromN(params.testLike,N_all);
%             end
%             X_test = data(test_indices,:);          
%             ll = likefunc(X_test,choles,ge.obsVar);
%             test_like = [test_like;ll];
%             loadSamples = true;
%         end                
% 
%     else        
%         load(params.initCond)
%         if isfield(state_sequence{end},'estimated_cholesky')
%             choles = state_sequence{end}.estimated_cholesky;
%             cc = cholcell(choles);
%         else
%             cc = state_sequence{end}.estimated_components;
%             choles = cellchol(cc,'cheat',true);
%         end
%         ge.cc = cc;   
%         X_test = data(test_indices,:);
%     end                   
%     
%     if params.sigmaSteps > 0 || params.startWithSigma
%         sigset_indices = chooseKfromN(batchSize,N_all);
%         X_sigset = data(sigset_indices,:);
%     end
%     
%     if params.startWithSigma
%         ge.obsVar = gestaltFindSigmaX(ge,choles,X_sigset,params.priorSamples,like_method,loadSamples,verb);
%         loadSamples = true;
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     for batch = 1:batchNum
%         if params.verbose == 1
%             printCounter(batch,'maxVal',batchNum,'stringVal','Batch');
%         elseif params.verbose >= 2
%             fprintf('Batch %d/%d ',batchNum,batch);
%         end
%         % sample a data batch
%         img_indices = chooseKfromN(batchSize,N_all);
%         X = data(img_indices,:);        
%         act_batch_like = [];
%         if ~strcmp(params.likeComp,'none')
%             ll = likefunc(X,choles,ge.obsVar);
%             act_batch_like = [act_batch_like ll];
%             loadSamples = true;
%         end        
%         
%         for step = 1:stepNum
%             if params.verbose == 2
%                 printCounter(step,'maxVal',stepNum,'stringVal','Step');
%             elseif params.verbose >= 3
%                 fprintf('Step %d/%d ',stepNum,step);
%             end
%             % calculate gradient
%             grad = gestaltLogLikelihoodGradient(ge,params.priorSamples,X,choles,'loadSamples',loadSamples,'method','scinot','template',t,'verbose',verb);            
%             % use the same set of samples at every iteration
%             loadSamples = true;
%             % update params 
%             choles = celladd(choles,1,grad,params.learningRate);                                    
%             
%             % calculate likelihood on batch
%             if ~strcmp(params.likeComp,'none')
%                 ll = likefunc(X,choles,ge.obsVar);
%                 act_batch_like = [act_batch_like ll];
%             end
%         end
%         batch_like = [batch_like; act_batch_like];
%         if strcmp(params.likeComp,'full')
%             ll = likefunc(data,choles,ge.obsVar);
%             full_like = [full_like;ll];
%         end
%         if ~isempty(test_indices)       
%             ll = likefunc(X_test,choles,ge.obsVar);
%             test_like = [test_like;ll];
%         end
%         % save stuff
%         state.estimated_components = cholcell(choles);
%         state.estimated_cholesky = choles;
%         state.sigma_x = ge.obsVar;
%         state.learningRate = params.learningRate;
%         state_sequence{end+1} = state;
%         batch_indices = [batch_indices; img_indices];
%         save('bin/gradasc_iter.mat','batch_like','full_like','test_like','ge','state_sequence', ...
%             'batch_indices','test_indices','sigset_indices','trueCC','trueSigma','-v7.3');
%         
%         if params.sigmaSteps > 0 && rem(batch,params.sigmaSteps) == 0
%             ge.obsVar = gestaltFindSigmaX(ge,choles,X_sigset,params.priorSamples,like_method,loadSamples,verb);
%         end
%     end    
% end