function ge = gestaltIEM(ge,X,nSamples,maxStep,lrate)
    ccInit = randomCovariances(ge.k,ge.Dv);
    X_old = ge.X;
    N_old = ge.N;
    cholesky = ccInit;
    for j=1:ge.k
        cholesky{j} = chol(cholesky{j});
    end
    ge.X = X;
    ge.N = size(ge.X,1);
    sdim = ge.k+(ge.Dv*ge.B);
    
    pCC{1} = ccInit;
    S = {};
    
    cc_next = cell(1,ge.k);
    samples = zeros(ge.N,nSamples,sdim);
    for i=1:maxStep
        cc_prev = ge.cc;
        for n=1:ge.N
            samples(n,:,:) = gestaltGibbs(ge,n,nSamples,'slice',0.05);            
            grad = gestaltParamGrad(ge,samples(n,:,:),cholesky);
            for j=1:ge.k
                cholesky{j} = cholesky{j} + lrate * grad{j};
                cc_next{j} = cholesky{j}' * cholesky{j};
            end
            ge.cc = cc_next;
        end
        
        diff = 0;
        for j=1:ge.k
            diff = diff + sum(sum((cc_next{j}-cc_prev{j})*(cc_next{j}-cc_prev{j})));           
        end
        
        pCC{i+1} = ge.cc;
        S{i} = samples;
        save('iter.mat','pCC','S');
    end
        
    ge.X = X_old;
    dnum = ge.N;
    ge.N = N_old;
    plotCovariances(ge,dnum,false);
end
