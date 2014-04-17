function ge = gestaltIEM(ge,X,nSamples,maxStep,lrate,precision,randseed)
    
    if strcmp(randseed,'last')
        load lastrandseed;
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');
            
    ccInit = randomCovariances(ge.k,ge.Dv,'precision',precision);
    
    X_old = ge.X;
    N_old = ge.N;
    
    if ~precision
        cc_old = ge.cc;
        ge.cc = ccInit;
    else
        cc_old = ge.pc;
        ge.pc = ccInit;
    end
    
    ge.X = X;
    ge.N = size(ge.X,1);
    sdim = ge.k+(ge.Dv*ge.B);
    % maximum change of a parameter over a cycle should not be more than:
    goaldiff = (1 / ge.N) * ones(ge.Dv);
    % empirical correction of the dimension dependence of the largest eigenvalue of the inverse covariance
    goaldiff = goaldiff / (ge.Dv * 0.025);
    
    cholesky = ccInit;
    for j=1:ge.k
        cholesky{j} = chol(cholesky{j});
    end    
    cholparnum = (ge.Dv^2 + ge.Dv) / 2;
    
    pCC{1} = ccInit;
    S = {};
    
    cc_next = cell(1,ge.k);
    samples = zeros(ge.N,nSamples,sdim);
    for i=1:maxStep
        fprintf('IEM cycle %d datapoint %d/',i,ge.N);
        if ~precision
            cc_prev = ge.cc;
        else
            cc_prev = ge.pc;
        end
        
        skipped = 0;
        avgrate = 0;
        for n=1:ge.N
            printCounter(n);
            fprintf(' ');
            
            % E-step: Gibbs sampling
            [samples(n,:,:),rr] = gestaltGibbs(ge,n,nSamples,0.05,'verbose',1,'precision',precision);            
            if rr < 0                
                fprintf('\b');                
                skipped = skipped + 1;
                continue;
            end
            
            % M-step: gradient ascent            
            grad = gestaltParamGrad(ge,samples(n,:,:),cholesky,'precision',precision);                        
            
            % choose learning rate
            meanvals = zeros(1,j);
            for j=1:ge.k
                meanvals(1,j) = meanvals(1,j) + mean(mean(abs(grad{j}),2),1);
            end
            meanval = mean(meanvals,2);   

            for j=1:ge.k
                % choose learning rate
                %actrate = min(goaldiff ./ abs(grad{j}),lrate * ones(ge.Dv));
                %actrate = min(goaldiff / meanval,lrate * ones(ge.Dv));
                actrate = min(goaldiff / meanvals(1,j),lrate * ones(ge.Dv));
                avgrate = avgrate + sum(sum(actrate))/cholparnum;
                % update 
                cholesky{j} = cholesky{j} + actrate .* grad{j};
                cc_next{j} = cholesky{j}' * cholesky{j};
            end
            
            % update parameters
            if ~precision
                ge.cc = cc_next;
            else
                ge.pc = cc_next;
            end
            
            for b=1:9+2*(floor(log10(nSamples))+1)
                fprintf('\b');
            end
        end
        
        diff = 0;
        for j=1:ge.k
            diff = diff + sum(sum((cc_next{j}-cc_prev{j})*(cc_next{j}-cc_prev{j})));           
        end
        diff = diff / (ge.k*ge.Dv^2);
        
        if ~precision
            pCC{i+1} = ge.cc;
        else
            pCC{i+1} = ge.pc;
        end
        
        S{i} = samples;
        save('iter.mat','pCC','S');
        fprintf(' avglr %.2e diff %.2e skipped %d\n',avgrate/(ge.N*ge.k),diff,skipped);
    end
        
    ge.X = X_old;
    dnum = ge.N;
    ge.N = N_old;
    
    if ~precision
        ge.pCC = ge.cc;        
        ge.cc = cc_old;
    else
        ge.pPC = ge.pc;        
        ge.pc = cc_old;
    end
    plotCovariances(ge,dnum,precision);
end
