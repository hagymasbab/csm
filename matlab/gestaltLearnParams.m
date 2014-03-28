function ge = gestaltLearnParams(ge,ccInit,X,nSamples,maxStep,varargin)
    p = inputParser;
    addParamValue(p,'test',false,@islogical);
    addParamValue(p,'close',false,@islogical);
    addParamValue(p,'firstSamples',{});
    p.KeepUnmatched = true;
    parse(p,varargin{:});
    test = p.Results.test;
    close = p.Results.close;   
    fSamp = p.Results.firstSamples;
    numGivenSamp = size(fSamp,2);

    if close
        fprintf('Test mode: using real covariance components plus noise as initial condition instead of the specified values!\n');
        ccInit = ge.cc;
        rancc = randomCovariances(ge.k,ge.Dv);
        for i=1:ge.k
            ccInit{i} = ccInit{i} + rancc{i};
        end
    end
    
    % store existing specifications from ge to restore them at the end
    cc_old = ge.cc;
    X_old = ge.X;
    N_old = ge.N;
    
    % put specified values to the model structure 
    ge.cc = ccInit;
    ge.X = X;
    ge.N = size(ge.X,1);
    
    % structures for saving iterations
    pCC{1} = ccInit;
    S = {};
    gVV = {};
    
    for i=1:maxStep
        fprintf('EM iteration %d\n..E-step\n....',i);
        % E-step: sampling the posterior
        if test
            fprintf('test mode: using data instead of samples!\n');
            samples = ones(ge.N,nSamples,ge.k+ge.Dv);
            for l=1:nSamples
                samples(:,l,1:ge.k) = ge.G;
                samples(:,l,ge.k+1:ge.k+ge.Dv) = ge.V;
            end
        elseif i<=numGivenSamp
            fprintf('using specified sample set\n');
            samples = fSamp{i};
        else
            samples = gestaltSamplePosterior(ge,nSamples);            
        end 
        
        % M-step: updating the parameters
        fprintf('..M-step\n');       
        VV = zeros(ge.Dv,ge.Dv);        

        fprintf('....Sample %d/', ge.N*nSamples);
        scalars = zeros(1,ge.k);
        matrices = zeros(ge.k,ge.Dv,ge.Dv);
        gVV{i} = cell(1,ge.k);
        for j=1:ge.k
            gVV{i}{j} = zeros(ge.Dv,ge.Dv);
        end
        for n=1:ge.N
            for l=1:nSamples
                printCounter((n-1)*nSamples+l);
                g = reshape(samples(n,l,1:ge.k),1,ge.k);
                v = reshape(samples(n,l,ge.k+1:ge.k+ge.Dv),1,ge.Dv);
                scalars = scalars + g .* g;
                vv = v'*v;
                VV = VV + vv;
                
                for j=1:ge.k
                    cc_temp = ge.cc;
                    cc_temp{j} = zeros(ge.Dv,ge.Dv);
                    gVV{i}{j} = gVV{i}{j} + g(1,j) * vv;
                    actmat = g(1,j) * (vv - componentSum(g',cc_temp));
                    matrices(j,:,:) = matrices(j,:,:) + reshape(actmat,1,ge.Dv,ge.Dv);
                end
            end
        end
        cc_next = cell(1,ge.k);
        for j=1:ge.k
            cc_next{j} = (1/scalars(1,j)) * reshape(matrices(j,:,:),ge.Dv,ge.Dv);
            % TEST - cheating !!! making the matrix pos def 
            L = ldl(cc_next{j});
            cc_next{j} = L*L';
            gVV{i}{j} = (1/scalars(1,j)) * gVV{i}{j};
        end        
        
        fprintf('\n');
        diff = 0;
        for j=1:ge.k        
            diff = diff + sum(sum((cc_next{j}-ge.cc{j})*(cc_next{j}-ge.cc{j})));
        end  
        diff = diff / (ge.k*ge.Dv^2);
        
        ge.cc = cc_next;        
        fprintf('..Average of squared differences of parameters between last two iterations: %e\n',diff);
        
        % save some data to disk
        pCC{i+1} = ge.cc;
        if ~test
            S{i} = samples;
        end
        VC{i} = VV / (nSamples*ge.N);
        save('iter.mat','pCC','S','VC','gVV');
        
        % test convergence
        if diff < 1e-3
            fprintf('Convergence achieved in %d steps.\n',i);
            break;
        end
    end
    
    % store results
    ge.pCC = ge.cc;
    
    % restore previously set fields in the model structure
    ge.cc = cc_old;
    ge.X = X_old;
    dnum = ge.N;
    ge.N = N_old;
    % visualise convergence    
    plotCovariances(ge,dnum);
end
                    