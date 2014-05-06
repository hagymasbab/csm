function [diff,longdiff] = gestaltIEM(ge,X,nSamples,maxStep,randseed,varargin)
    parser = inputParser;
    addParamValue(parser,'learningRate',0.01,@isnumeric);
    addParamValue(parser,'rateMethod','componentwise_goal');
    addParamValue(parser,'plot',1,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    addParamValue(parser,'multistep',false,@islogical);
    addParamValue(parser,'verbose',2,@isnumeric);
    addParamValue(parser,'approximatePostCov',false,@islogical);
    parse(parser,varargin{:});
    lrate = parser.Results.learningRate; 
    ratemethod = parser.Results.rateMethod; 
    plot = parser.Results.plot;
    approx = parser.Results.approximatePostCov;
    precision = parser.Results.precision;
    verb = parser.Results.verbose;
    multistep = parser.Results.multistep;
    
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
    
    minperm = [];
    diff = zeros(1,maxStep+1);
    diff(1,1) = covcompRootMeanSquare(ccInit,cc_old,minperm);
    if nargout > 1
        longdiff = zeros(1,maxStep*ge.N+1);
        longdiff(1,1) = diff(1,1);
    end    
    
    ge.X = X;
    ge.N = size(ge.X,1);
    sdim = ge.k+(ge.Dv*ge.B);
    % maximum change of a parameter over a cycle should not be more than:
    goaldiff = (1 / ge.N) * ones(ge.Dv);
    % empirical correction of the dimension dependence of the largest eigenvalue of the inverse covariance
    goaldiff = goaldiff / (ge.Dv * 0.025);
    
    cholesky = cellchol(ccInit);  
    cholparnum = (ge.Dv^2 + ge.Dv) / 2;
    
    pCC{1} = ccInit;
    S = {};
    
    if plot>1
        subplot = @(m,n,p) subtightplot (m, n, p, [0.025 0.001], [0 0.025], [0 0.01]);
        clf;
    end
    
    cc_next = cell(1,ge.k);
    samples = zeros(ge.N,nSamples,sdim);
    for i=1:maxStep
        if verb > 1
            fprintf('IEM cycle %d datapoint %d/',i,ge.N);
        end
        if plot>1
            nopause = false;
        end
        if ~precision
            cc_prev = ge.cc;
        else
            cc_prev = ge.pc;
        end
        
        skipped = 0;
        avgrate = 0;
        for n=1:ge.N
            if verb>1
                printCounter(n);
            end
            if verb > 1
                fprintf(' ');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % ACTUAL IEM PART BEGINS HERE            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % E-step: Gibbs sampling
            [samples(n,:,:),rr] = gestaltGibbs(ge,n,nSamples,'verbose',verb-1,'precision',precision,'approximatePostCov',approx);            
            if rr < 0                
                if verb>1
                    fprintf('\b');                
                end
                skipped = skipped + 1;
                continue;
            end
            
            % M-step: gradient ascent            
            grad = gestaltParamGrad(ge,samples(n,:,:),cholesky,'precision',precision);                        
            
            % choose learning rate
            meanvals = zeros(1,ge.k);
            for j=1:ge.k
                meanvals(1,j) = meanvals(1,j) + mean(mean(abs(grad{j}),2),1);
            end
            meanval = mean(meanvals,2);   
            
            chol_cand = cholesky;
            cdll = gestaltCompleteDataLogLikelihood(ge,samples(n,:,:),cholesky);
            actrate = cell(1,ge.k);
            avgrates = zeros(1,ge.k);
            for j=1:ge.k
                % choose learning rate
                if strcmp(ratemethod,'elementwise_goal')
                    actrate = min(goaldiff ./ abs(grad{j}),lrate * ones(ge.Dv));
                elseif strcmp(ratemethod,'mean_goal')
                    actrate = min(goaldiff / meanval,lrate * ones(ge.Dv));
                elseif strcmp(ratemethod,'componentwise_goal')                  
                    actrate{j} = min(goaldiff / meanvals(1,j),lrate * ones(ge.Dv));
                elseif strcmp(ratemethod,'only_goal')
                    actrate{j} = goaldiff / meanvals(1,j);
                elseif strcmp(ratemethod,'only_rate')
                    actrate{j} = lrate * ones(ge.Dv);
                else
                    exit('Invalid learning rate method: %s',ratemethod);
                end
                avgrates(1,j) = sum(sum(actrate{j}))/cholparnum;
            
                % update 
                if ~multistep
                    cholesky{j} = cholesky{j} + actrate{j} .* triu(grad{j});
                end
            end    
            if multistep
                increase = true;
                while increase
                    for j=1:ge.k
                        chol_cand{j} = chol_cand{j} + actrate{j} .* triu(grad{j});
                    end
                    cdll_cand = gestaltCompleteDataLogLikelihood(ge,samples(n,:,:),chol_cand);
%                     fprintf(' %f %f\n',cdll,cdll_cand);
%                     pause;
%                     for b=1:3+2*7
%                         fprintf('\b');
%                     end
                    if cdll_cand > cdll
                        cholesky = chol_cand;
                        cdll = cdll_cand;
                    else
                        increase = false;
                    end
                end
            end
            for j=1:ge.k
                cc_next{j} = cholesky{j}' * cholesky{j};                                             
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % ACTUAL IEM PART ENDS HERE            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if plot>1
                [nopause,plot] = plotIEMStep(ge,cholesky,i,n,samples,grad,cc_next,cc_old,avgrates,nopause);
            end
            
            if nargout > 1
                lidx = 1+(i-1)*ge.N+n;
                [longdiff(1,lidx),minperm] = covcompRootMeanSquare(cc_next,ge.cc,minperm);
            end
            
            % update parameters
            if ~precision
                ge.cc = cc_next;
            else
                ge.pc = cc_next;
            end
            
            if verb>1
                for b=1:9+2*(floor(log10(nSamples))+1)
                    fprintf('\b');
                end
            end
        end
        
        reldiff = covcompRootMeanSquare(cc_next,cc_prev,1:ge.k);
        [diff(1,i+1),minperm] = covcompRootMeanSquare(cc_next,cc_old,minperm);        
        
        if ~precision
            pCC{i+1} = ge.cc;
        else
            pCC{i+1} = ge.pc;
        end
        
        S{i} = samples;
        save('iter.mat','pCC','S');
        if verb>1
            fprintf(' avglr %.2e diff %.2e skipped %d\n',sum(avgrates,2)/(ge.N*ge.k),reldiff,skipped);
        end
        
        if reldiff < 1e-3
            if verb>1
                fprintf('Convergence achieved in %d steps.\n',i);
            end
            break;
        end
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
    if plot>0
        plotCovariances(ge,dnum,precision);
    end
end
    
function [nopause,plot] = plotIEMStep(ge,cholesky,i,n,samples,grad,cc_next,cc_old,avgrates,nopause)
    hor = 6;
    for j=1:ge.k
        % cholesky
        subplot(ge.k,hor,(j-1)*hor+1);
        viewImage(cholesky{j},'magnif',false);
        title(sprintf('chol %d at %d#%d',j,i,n));                    
        % gradients
        subplot(ge.k,hor,(j-1)*hor+3);
        %compgrad = grad{j}/oldchol{j};
        compgrad = grad{j};
        viewImage(compgrad,'magnif',false);
        gavg = mean(squeeze(samples(n,:,j)));
        title(sprintf('grad sg%d=%.3f',j,gavg));
        % delta
        subplot(ge.k,hor,(j-1)*hor+4);
        viewImage(compgrad*avgrates(1,j),'magnif',false);                    
        title(sprintf('delta lr=%.3f',avgrates(1,j)));
        % next components
        subplot(ge.k,hor,(j-1)*hor+5);
        viewImage(cc_next{j},'magnif',false);
        title(sprintf('comp %d at %d#%d',j,i,n));
        % truth
        subplot(ge.k,hor,(j-1)*hor+6);
        viewImage(cc_old{j},'magnif',false);
        title(sprintf('true comp %d',j));
    end
    % data cov                    
    subplot(ge.k,hor,2);
    viewImage(cov(squeeze(ge.X(n,:,:))),'magnif',false);
    title(sprintf('data cov, g1=%.3f',ge.G(n,1)));
    % sample cov                    
    subplot(ge.k,hor,hor+2);
    vsamp = reshape(samples(n,:,ge.k+1:ge.k+(ge.Dv*ge.B)),size(samples,2)*ge.B,ge.Dv);
    viewImage(cov(vsamp),'magnif',false);
    title('sample cov');
    pause(0.01);
    plot = 2;    
    if ~nopause
        ch = getkey('non-ascii');
        if strcmp('f',ch)
            nopause = true;
        elseif strcmp('r',ch)
            plot = 1;
        end
    end
end