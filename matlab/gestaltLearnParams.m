function ge = gestaltLearnParams(ge,init,X,nSamples,maxStep,varargin)
    p = inputParser;
    addParamValue(p,'test',false,@isnumeric);
    addParamValue(p,'learningRate',0.001,@isnumeric);
    addParamValue(p,'close',false,@islogical);
    addParamValue(p,'precision',false,@islogical);
    addParamValue(p,'firstSamples',{});
    p.KeepUnmatched = true;
    parse(p,varargin{:});
    test = p.Results.test;
    close = p.Results.close;   
    precision = p.Results.precision;   
    % scale learning rate with number of datapoints
    lrate = p.Results.learningRate / ge.N;
    fSamp = p.Results.firstSamples;
    numGivenSamp = size(fSamp,2);

    if strcmp(init,'close')
        % TODO implement precision component case
        fprintf('Using real components plus noise as initial condition.\n');
        ccInit = ge.cc;
        rancc = randomCovariances(ge.k,ge.Dv);
        for i=1:ge.k
            ccInit{i} = ccInit{i} + 5 * rancc{i};
        end
    elseif strcmp(init,'random')
        fprintf('Using random initial condition.\n');
        ccInit = randomCovariances(ge.k,ge.Dv,precision);
    end
    
    % store existing specifications from ge to restore them at the end
    X_old = ge.X;
    N_old = ge.N;
    
    if ~precision
        cc_old = ge.cc;    
        ge.cc = ccInit;
    else
        cc_old = ge.pc;    
        ge.pc = ccInit;   
    end
    
    cholesky = ccInit;
    for j=1:ge.k
        cholesky{j} = chol(cholesky{j});
    end
    ge.X = X;
    ge.N = size(ge.X,1);
    sdim = ge.k+(ge.Dv*ge.B);
    
    % structures for saving iterations
    pCC{1} = ccInit;
    S = {};
    gVV = {};
    cc_next = cell(1,ge.k);
    
    for i=1:maxStep
        fprintf('EM iteration %d: E-step - ',i);
        % E-step: sampling the posterior
        if test == 3
            fprintf('test mode: using data instead of samples!\n');
            samples = ones(ge.N,nSamples,sdim);
            randCov = randomCovariances(1,sdim);            
            noiseVar = (maxStep-i) * 1;
            randCov = noiseVar * randCov{1};
            for l=1:nSamples
                samples(:,l,1:ge.k) = ge.G(1:ge.N,:);
                samples(:,l,ge.k+1:sdim) = reshape(ge.V(1:ge.N,:,:),ge.N,ge.Dv*ge.B);
                samples(:,l,:) = mvnrnd(squeeze(samples(:,l,:)),randCov);
            end
        elseif i<=numGivenSamp
            fprintf('using specified sample set\n');
            samples = fSamp{i};
        else
            samples = gestaltSamplePosterior(ge,nSamples,'precision',precision);            
        end 
        
        if test == 1
            fprintf('....test mode: using ground truth for g instead of samples!\n');
            for l=1:nSamples
                samples(:,l,1:ge.k) = ge.G(1:ge.N,:);                
            end
        elseif test == 2
            fprintf('....test mode: using ground truth for v instead of samples!\n');
            for l=1:nSamples
                samples(:,l,ge.k+1:sdim) = reshape(ge.V(1:ge.N,:,:),1,ge.Dv*ge.B);
            end
        end
        
        % M-step: updating the parameters
        fprintf('\b M-step - ');       
        
        if ~precision
            grad = gestaltParamGrad(ge,samples,cholesky);
        else
            grad = gestaltParamGradPrec(ge,samples,cholesky);
        end
        
        % set the learning rate to change the matrices at most by about 0.1
        maxval = 0;
        for j=1:ge.k
            actmax = max(max(abs(grad{j})));
            if actmax > maxval
                maxval = actmax;
            end
        end
        lrate = min(0.15/maxval,lrate);
        if precision
            lrate = 1/lrate
        end
        
        for j=1:ge.k
            cholesky{j} = cholesky{j} + lrate * grad{j};
            cc_next{j} = cholesky{j}' * cholesky{j};
        end
        
        diff = 0;
        for j=1:ge.k
            if ~precision
                diff = diff + sum(sum((cc_next{j}-ge.cc{j})*(cc_next{j}-ge.cc{j})));
            else
                diff = diff + sum(sum((cc_next{j}-ge.pc{j})*(cc_next{j}-ge.pc{j})));
            end
        end  
        diff = diff / (ge.k*ge.Dv^2);
        
        if ~precision
            ge.cc = cc_next;        
        else
            ge.pc = cc_next;
        end
        fprintf('\b lr %.2e diff %.2e\n',lrate,diff);
        
        % save some data to disk
        pCC{i+1} = cc_next;
        if ~test
            S{i} = samples;
        end
        save('iter.mat','pCC','S');
        
        % test convergence
        if diff < 1e-3
            fprintf('Convergence achieved in %d steps.\n',i);
            break;
        end
    end
    
    % store results
    if ~precision    
        ge.pCC = ge.cc;        
        ge.cc = cc_old;
    else
        ge.pPC = ge.pc;        
        ge.pc = cc_old;
    end
    
    % restore previously set fields in the model structure
    ge.X = X_old;
    dnum = ge.N;
    ge.N = N_old;
    % visualise convergence    
    plotCovariances(ge,dnum,precision);
end
                    