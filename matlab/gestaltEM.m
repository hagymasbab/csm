function [cholesky,cc_next] = gestaltEM(ge,X,emBatchSize,maxStep,nSamples,randseed,varargin)        
    parser = inputParser;
    addParameter(parser,'learningRate',0,@isnumeric);    
    addParameter(parser,'plot',0,@isnumeric);
    addParameter(parser,'precision',false,@islogical);      
    addParameter(parser,'verbose',1,@isnumeric);
    addParameter(parser,'initCond','empty');
    addParameter(parser,'priorG','gamma');
    addParameter(parser,'sampler','gibbs');
    addParameter(parser,'syntheticData',true,@islogical);    
    addParameter(parser,'burnin',0,@isnumeric);
    addParameter(parser,'stoppingDiff',0,@isnumeric);      
    addParameter(parser,'computeLikelihood','none');      
    addParameter(parser,'likelihoodSamples',50,@isnumeric);   
    addParameter(parser,'cctComponents',false,@islogical); 
    addParameter(parser,'savingCode',0,@isnumeric);   
    addParameter(parser,'skipCheck',false,@islogical);
    parse(parser,varargin{:});        
    params = parser.Results;      
    
    % VERBOSITY LEVELS
    % 0 - don't print anything
    % 1 - only print EM-step number progression
    % 2 - print progress within EM steps, cycling over data and components
    % 3 - print progress in sampling too
    
    
    if params.plot>1
        % redefine plot to make figures more dense
        subplot = @(m,n,p) subtightplot (m, n, p, [0.025 0.001], [0 0.025], [0 0.01]);
        clf;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    % SET RANDOM SEED          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(randseed,'last')
        load lastrandseed;
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    % CREATE INITAL PARAMETER MATRICES AND MODEL STRUCTURE     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
    if params.savingCode == 0
        savingCode = randi(1000000);
    else
        savingCode = params.savingCode;
    end
    ge_saved = ge;

    if ischar(params.initCond)
        if strcmp(params.initCond,'random') 
            ccInit = randomCovariances(ge.k,ge.Dv,'precision',params.precision);        
        elseif strcmp(params.initCond,'shifted') 
            ccInit = gestaltCovariances(ge.k,ge.A',floor(sqrt(ge.Dx)/4));  
        elseif strcmp(params.initCond,'empty') 
            ccInit = cell(1,ge.k);
            for kk = 1:ge.k
                if ge.nullComponent && kk == ge.k
                    varScale = 1;
                else
                    varScale = 0.001;
                end
                ccInit{kk} = varScale * eye(ge.Dv);
            end
        elseif strcmp(params.initCond,'pretrain')
            pre_seed = randi(intmax);
            data = reshape(X,size(X,1)*size(X,2),size(X,3));
            ge.X = repmat(data,floor(1000/size(data,1)),1);
            ccInit = gestaltPretrain(ge,1000,pre_seed,'plot',false);
            for kk = 1:ge.k
                ccInit{kk} = 10 * ccInit{kk};
                ccInit{kk} = max(ccInit{kk}(:)) * eye(ge.Dv);
            end
        end
    elseif iscell(params.initCond) && size(params.initCond) == ge.k
        ccInit = params.initCond;
    elseif isnumeric(params.initCond) && length(params.initCond) == 1
        ccInit = cell(1,ge.k);
        savingCode = params.initCond;
        load(sprintf('cc_%d.mat',savingCode))
        for i=1:ge.k
            [~,e] = chol(cc_saved{i});
            if e == 0
                ccInit{i} = cc_saved{i};
            else
                ccInit{i} = nearestSPD(cc_saved{i});
            end
        end
    else
        fprintf('not implemented');
        return;
    end
    
    %fprintf('The saving code for this run is %d\n',savingCode);
    savename = sprintf('bin/cc_%d.mat',savingCode);
    
    cholesky = cellchol(ccInit);                  
    goal_cc = extractComponents(ge,params.precision);    
    ge = replaceComponents(ge,ccInit,params.precision);

    sdim = ge.k+(ge.Dv*ge.B);
    % correct learning rate for data set size
    %params.learningRate = params.learningRate * (10 / emBatchSize);
%     ge.N = size(ge.X,1);  
    if ndims(X) == 2
        if ge.B > 1
            throw(MException('Gestalt:EM:WrongBatchSize','The data is not organised to batches while the model is specified for batches of observarions of size %d',ge.B));
        end
        X = reshape(X,size(X,1),1,ge.Dv);
    end
    
    adjustLR = params.learningRate == 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    % INITIALISE ARRAYS FOR SAVING VALUES   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            
    
    truthdiff_diagonals = false;
    % the best permutation for comparing estimated components to truth    
    state_sequence = cell(1,maxStep+1);   
    true_c = componentSum(ones(ge.k,1),goal_cc);
    act_c = componentSum(ones(ge.k,1),ccInit);
    if params.syntheticData
        true_c = cov(squeeze(mean(ge.V,2))); % this is assuming that when we have synthetic data, we generated it using ge.V
        state.difference_to_truth = covcompRootMeanSquare(act_c,true_c,1,'useDiagonals',truthdiff_diagonals);        
%         if params.detailedDiff
%             state.detailed_difference = covcompRootMeanSquare(ccInit,goal_cc,[],'verbose',true);
%         end
    end
    state.matrix_norms = {};
    for i=1:ge.k
        state.matrix_norms{i} = norm(ccInit{i});
    end
    state.relative_difference = 0;        
    state.estimated_components = ccInit;
    state.samples = [];    
    
    if strcmp(params.computeLikelihood,'full')
        % likelihoods for real images tend to be smaller, so we have to
        % calculate the likelihood with recording coefficients and
        % exponents at every step in that case
        % use_scinot = ~params.syntheticData;
        % calculate log-likelihood on full dataset
        state.full_loglike = gestaltLogLikelihood(ge,params.likelihoodSamples,X);
        best_loglike = state.full_loglike;
        if params.verbose>1
            fprintf('Log-likelihood on full dataset: %f\n',state.full_loglike);
        end
    end
    
    state_sequence{1} = state;                  
            
    cc_next = cell(1,ge.k);
    
    for step=1:maxStep
        if params.verbose == 1
            printCounter(step,'maxVal',maxStep,'stringVal','EM step')
        end
        if params.verbose > 1
            fprintf('EM step %d/%d, saving code %d',maxStep,step,savingCode)
        end
        if params.verbose == 2
            printProgress(emBatchSize,'Datapoint');
        end
        
        cc_prev = extractComponents(ge,params.precision);
             
        % choose a batch from tha big dataset randomly
        img_indices = chooseKfromN(emBatchSize,size(X,1));
        ge.X = X(img_indices,:,:);        
        
        % check whether current components are sparse         
        ge.sparseComponents = length(find(componentSum(1,ge.cc))) < ge.Dv^2 / 2;        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % E - step            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        vsamp = zeros(emBatchSize,nSamples,ge.B,ge.Dv);
        gsamp = zeros(emBatchSize,nSamples,ge.k);
        skipped = 0;
        parfor n=1:emBatchSize
        %for n=1:emBatchSize
            
             if params.verbose == 3 && ~strcmp(params.sampler,'test')
                 fprintf('\nDatapoint %d/%d ',emBatchSize,n);                   
             elseif params.verbose ==2    
                 printProgress()
             end

            % sampling
            initG = (1/ge.k) * ones(ge.k,1);
            try
                if strcmp(params.sampler,'gibbs')
                    [vsamp(n,:,:,:),gsamp(n,:,:),act_zsamp,~] = gestaltGibbs(ge,n,nSamples,'verbose',params.verbose-2,'precision',params.precision, ...
                        'initG',initG,'contrast',ge.contrast,'burnin',params.burnin,'skipCheck',params.skipCheck);            
                    %mean(gsamp(n,:,:))
                    %mean(act_zsamp)
                elseif strcmp(params.sampler,'test')
                     % this is assuming that when we have synthetic data,
                     % we generated it using ge.V, ge.G and ge.Z
                     vsamp(n,:,:,:) = permute(repmat(reshape(ge.V(n,:,:),ge.B,ge.Dv),[1,1,nSamples]),[3,1,2]);
                     gsamp(n,:,:) = repmat(ge.G(n,:),[nSamples,1]);
                end                    
            catch e
                % if couldn't find a valid g-sample in 10 steps, skip                
                if (strcmp(e.identifier,'Gestalt:Gibbs:TooManyTries'))
                    skipped = skipped + 1;
                    continue;
                else
                    rethrow(e);
                end
            end
        end        
        
        if strcmp(params.computeLikelihood,'batch')
            state.batch_loglike_pre = gestaltLogLikelihood(ge,params.likelihoodSamples,ge.X);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % M - step            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % gradient of the parameters of the complete-data log-likelihood            
        grad = gestaltParamGrad(ge,vsamp,gsamp,cholesky,'precision',params.precision,'verbose',params.verbose-1);                        
                
        if adjustLR
%             products = cell(1,ge.k);
%             for kk = 1:ge.k
%                 products{kk} = grad{kk}*grad{kk}';
%             end            
%             matgrad = cell2mat(products);
%             maxsq = sqrt(max(matgrad(:)));
%             
            matgrad = cell2mat(grad);
            maxsq = max(matgrad(:));
            
            max_LR = 10^(1-floor(log10(maxsq)+log10(ge.Dv)));                        
            if step == 1
                params.learningRate = max_LR;
            else
%                 params.learningRate = min(params.learningRate,max_LR);
            end
        end
        
        %params.learningRate
        
        % update cholesky components
        for j=1:ge.k
            if ge.nullComponent && j==ge.k
                cc_next{j} = ccInit{j};
            else
                cholesky{j} = cholesky{j} + params.learningRate .* grad{j};
                cc_next{j} = cholesky{j}' * cholesky{j};                                             
            end
        end

        ge = replaceComponents(ge,cc_next,params.precision);      
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % PRINT AND SAVE DATA            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
                
        [state.relative_difference,~,maxel_diff_rel] = covcompRootMeanSquare(cc_next,cc_prev,1:ge.k);
        act_c = componentSum(ones(ge.k,1),cc_next);
        if params.syntheticData
            [state.difference_to_truth,~,maxel_diff] = covcompRootMeanSquare(act_c,true_c,1,'useDiagonals',truthdiff_diagonals);
        end
        state.estimated_components = extractComponents(ge,params.precision);
        state.img_indices = img_indices;
        %state.samples = samples;
        state.matrix_norms = {};
        for i=1:ge.k
            state.matrix_norms{i} = norm(cc_next{i});
        end                    
        state.learningRate = params.learningRate;
                
        if params.verbose > 1
            synstr = '';
            if params.syntheticData
                synstr = sprintf(' truth_offd %f maxtruth_offd %f',state.difference_to_truth,maxel_diff);
            end
            fprintf(' maxreldiff %.2e %s skipped %d\n',maxel_diff_rel,synstr,skipped);
        end
        
        if strcmp(params.computeLikelihood,'full')
            % calculate log-likelihood on full dataset
            state.full_loglike = gestaltLogLikelihood(ge,params.likelihoodSamples,X);
            best_so_far = false;
            if state.full_loglike > best_loglike
                best_loglike = state.full_loglike;
                best_so_far = true;
            end
            if params.verbose>1
                fprintf('Log-likelihood on full dataset: %f\n',state.full_loglike);
            end
        elseif strcmp(params.computeLikelihood,'batch')
            state.batch_loglike = gestaltLogLikelihood(ge,params.likelihoodSamples,ge.X);
            if params.verbose>1
                fprintf('Log-likelihood on actual batch before learning: %f\n',state.batch_loglike_pre);
                fprintf('                               after learning: %f\n',state.batch_loglike);
            end
        end
        
        state_sequence{step+1} = state; 
        
        if params.syntheticData
            save('bin/iter.mat','state_sequence','goal_cc','ge_saved','-v7.3');
        else
            save('bin/iter.mat','state_sequence','ge_saved','-v7.3');            
        end
        
        if ~strcmp(params.computeLikelihood,'full') || best_so_far
            cc_saved = cc_next;        
            save(savename,'cc_saved','ge_saved');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % TEST FOR CONVERGENCE   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        % if the largest (including the diagonal) change is less than one over ten
        % thousand, we are safe to stop
        if params.stoppingDiff > 0 && maxel_diff_rel < params.stoppingDiff
            if params.verbose>1
                fprintf('Convergence achieved in %d steps.\n',step);
            end
            break;
        end
        
        %whos
    end        

    if params.plot>0
        if ge.k > 5
            plotsize = ceil(sqrt(ge.k));
            for act_k=1:ge.k
                subplot(plotsize,plotsize,act_k);
                viewImage(cc_next{act_k});
            end
        else
            ge = replaceComponents(ge,goal_cc,params.precision);
            plotCovariances(ge,ge.N,params.precision,[]);
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
% FUNCTIONS FOR CALCULATIONS, PLOTTING AND DATA ACCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ge = replaceComponents(ge,comps,precision)
    if ~precision       
        ge.cc = comps;
    else      
        ge.pc = comps;
    end
end

function comps = extractComponents(ge,precision)
    if ~precision
        comps = ge.cc;        
    else
        comps = ge.pc;        
    end
end
