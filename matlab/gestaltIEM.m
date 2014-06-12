function [diff,like] = gestaltIEM(ge,X,nSamples,maxStep,randseed,varargin)        
    parser = inputParser;
    addParamValue(parser,'learningRate',0.01,@isnumeric);
    addParamValue(parser,'rateMethod','componentwise_goal');
    addParamValue(parser,'plot',1,@isnumeric);
    addParamValue(parser,'precision',false,@islogical);
    addParamValue(parser,'multistep',false,@islogical);
    addParamValue(parser,'verbose',2,@isnumeric);
    addParamValue(parser,'calculateLikelihood',false,@islogical);
    addParamValue(parser,'increaseLikelihood',true,@islogical);
    addParamValue(parser,'likelihoodSamples',10,@isnumeric);
    addParamValue(parser,'fullLikelihood',true,@islogical);
    addParamValue(parser,'noiseLevel',0.1,@isnumeric);
    addParamValue(parser,'sortData',false,@islogical);
    parse(parser,varargin{:});
    
    lrate = parser.Results.learningRate; 
    ratemethod = parser.Results.rateMethod; 
    plot = parser.Results.plot;
    precision = parser.Results.precision;
    verb = parser.Results.verbose;
    multistep = parser.Results.multistep;
    calcLike = parser.Results.calculateLikelihood;
    incLike = parser.Results.increaseLikelihood;
    fullLike = parser.Results.fullLikelihood;
    likeSamp = parser.Results.likelihoodSamples;    
    noiseLevel = parser.Results.noiseLevel;   
    sortData = parser.Results.sortData;   
    
    if incLike || fullLike
        calcLike = true;
    end
    
    if plot>1
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
        
    ccInit = randomCovariances(ge.k,ge.Dv,'precision',precision,'noiseLevel',noiseLevel);        
    cholesky = cellchol(ccInit);   
    pCC{1} = ccInit;            
    
    cc_old = extractComponents(ge,precision);    
    ge = replaceComponents(ge,ccInit,precision);
    ge.X = X;
    ge.N = size(ge.X,1);
    if sortData
        % here's a strong assumption that we used the first part of the
        % data stored in the structure
        ge.G = ge.G(1:ge.N,:);
        ge = sortByGestalt(ge);
    end    
    sdim = ge.k+(ge.Dv*ge.B);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    % INITIALISE ARRAYS FOR SAVING VALUES   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
    
    minperm = [];
    diff = zeros(1,maxStep+1);
    diff(1,1) = covcompRootMeanSquare(ccInit,cc_old,minperm);
    longdiff = zeros(1,maxStep*ge.N+1);
    if nargout > 1        
        longdiff(1,1) = diff(1,1);
    end              
        
    loglike = zeros(1,maxStep*ge.N+1);    
    like = zeros(1,maxStep+1);    
    if calcLike
        like(1) = gestaltLogLikelihood(ge,likeSamp,0,[]);
        loglike(1) = like(1);
    end
    
    S = {};    
    samples = zeros(ge.N,nSamples,sdim);
    mean_gradient = zeros(1,maxStep*ge.N);
    cc_next = cell(1,ge.k);
    
    for i=1:maxStep
        if verb == 2
            fprintf('IEM cycle %d datapoint %d/',i,ge.N);            
        end
        if plot>1
            nopause = false;
        end
        
        cc_prev = extractComponents(ge,precision);
        skipped = 0;
        avgrate = 0;
        
        for n=1:ge.N
            lidx = 1+(i-1)*ge.N+n;
            if verb==2
                printCounter(n);
                fprintf(' ');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % E - step            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Gibbs sampling
            [samples(n,:,:),rr] = gestaltGibbs(ge,n,nSamples,'verbose',verb-1,'precision',precision);            
            % if couldn't find a valid g-sample in 10 steps, skip
            if rr < 0                %%%%
                if verb==2
                    if rr == -1
                        fprintf('\b');                
                    else
                        delPrint(-rr-1);
                    end
                end
                skipped = skipped + 1;
                continue;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % M - step            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % gradient of the parameters of the complete-data log-likelihood            
            grad = gestaltParamGrad(ge,samples(n,:,:),cholesky,'precision',precision);                        
            
            % calculate rates and average gradients
            [actrate,avgrates,gradAvg] = rateMatrices(grad,ratemethod,lrate,ge.N);            
            avgrate = avgrate + sum(avgrates,2)/ge.k;    
            mean_gradient(1,(i-1)*ge.N+n) = gradAvg;    
            
            % update cholesky components
            old_chol = cholesky;
            if ~multistep
                for j=1:ge.k
                    % TEST
                    if n==1
                        actrate{j} = 0.1 * ones(ge.Dv);
                    end
                    cholesky{j} = cholesky{j} + actrate{j} .* triu(grad{j});
                end
            else                        
                chol_cand = cholesky;
                %cdll = gestaltCompleteDataLogLikelihood(ge,samples(n,:,:),cholesky);            
                %cdll = gestaltLogLikelihood(ge,5,0,cholesky);            
                cdll = gestaltLogLikelihood(ge,5,n,cholesky);            
                increase = true;
                while increase
                    for j=1:ge.k
                        chol_cand{j} = chol_cand{j} + actrate{j} .* triu(grad{j});
                    end
                    %cdll_cand = gestaltCompleteDataLogLikelihood(ge,samples(n,:,:),chol_cand);
                    %cdll_cand = gestaltLogLikelihood(ge,5,0,chol_cand);
                    cdll_cand = gestaltLogLikelihood(ge,5,n,chol_cand);

                    if cdll_cand > cdll
                        cholesky = chol_cand;
                        cdll = cdll_cand;
                    else
                        increase = false;
                    end
                end
            end            
            % update component matrices
            for j=1:ge.k
                cc_next{j} = cholesky{j}' * cholesky{j};                                             
            end
            cc_temp = extractComponents(ge,precision);
            ge = replaceComponents(ge,cc_next,precision);      
            
            if ge.k == 1 && ge.Dv == 2 && verb == 3
                dcov = cov(squeeze(ge.X(n,:,:)));
                fprintf('Data cov %.2f %.2f %.2f Chol %.2f %.2f %.2f Grad %.2f %.2f %.2f\n',dcov(1,1),dcov(2,2),dcov(1,2),cholesky{1}(1,1),cholesky{1}(2,2),cholesky{1}(1,2),grad{1}(1,1),grad{1}(2,2),grad{1}(1,2))
                pause;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % CONTROL FOR LIKELIHOOD            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            
            if calcLike
                if fullLike
                    loglike(lidx) = gestaltLogLikelihood(ge,likeSamp,0,[]);
                else
                    loglike(lidx) = gestaltLogLikelihood(ge,likeSamp,n,[]);
                end
                if incLike
                    % if likelihood didn't increase, revert
                    %if loglike(lidx) < loglike(lidx-1)
                    % TEST
                    if loglike(lidx) < loglike(lidx-1) && n > 1
                        cholesky = old_chol;
                        ge = replaceComponents(ge,cc_temp,precision);
                        skipped = skipped + 1;
                        loglike(lidx) = loglike(lidx-1);         
                        if ge.k == 1 && ge.Dv == 2 && verb == 3
                            fprintf('\b skipped\n');
                        end
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % PLOT, PRINT AND SAVE DATA            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                   
            
            if plot>1
                [nopause,plot] = plotIEMStep(ge,cholesky,i,n,samples,grad,cc_next,cc_old,avgrates,nopause);
            end
            
            if nargout > 1
                [longdiff(1,lidx),minperm] = covcompRootMeanSquare(cc_next,cc_old,[]);
            end                        
            
            if verb==2
                delPrint(nSamples);
            end
        end %for n=1:ge.N
        avgrate = avgrate / ge.N;
        
        reldiff = covcompRootMeanSquare(cc_next,cc_prev,1:ge.k);
        [diff(1,i+1),minperm] = covcompRootMeanSquare(cc_next,cc_old,[]);   
        if fullLike
            like(1,i+1) = loglike(1,1+i*ge.N);
        else
            like(1,i+1) = gestaltLogLikelihood(ge,likeSamp,0,[]);
        end
        
        pCC{i+1} = extractComponents(ge,precision);
        
        S{i} = samples;
        save('iter.mat','pCC','S','mean_gradient','diff','longdiff','loglike','like');
        if verb == 2
            fprintf(' avglr %.2e diff %.2e skipped %d\n',avgrate,reldiff,skipped);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        % TEST FOR CONVERGENCE   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %reldiff = covcompRootMeanSquare(cc_next,cc_prev,1:ge.k);
        if reldiff < 1e-6
            if verb>1
                fprintf('Convergence achieved in %d steps.\n',i);
            end
            break;
        end
    end        

    if plot>0
        ge = replaceComponents(ge,cc_old,precision);
        plotCovariances(ge,ge.N,precision,[]);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
% FUNCTIONS FOR CALCULATIONS, PLOTTING AND DATA ACCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [actrate,avgrates,gradAvg] = rateMatrices(grad,ratemethod,lrate,N)
    Dv = size(grad{1},1);
    k = size(grad,2);

    % mean of gradient for learning rate calculation
    [gradAvg,gradCompAvgs] = componentAbsMeans(grad);     
    
    % maximum change of a parameter over a cycle should not be more than:
    goaldiff = (1 / N) * ones(Dv);
    % empirical correction of the dimension dependence of the largest eigenvalue of the inverse covariance
    goaldiff = goaldiff / (Dv * 0.025);
    
    actrate = cell(1,k);
    avgrates = zeros(1,k);
    for j=1:k
        % choose learning rate according to selected method
        if strcmp(ratemethod,'elementwise_goal')
            actrate{j} = min(goaldiff ./ abs(grad{j}),lrate * ones(Dv));
        elseif strcmp(ratemethod,'mean_goal')
            actrate{j} = min(goaldiff / gradAvg,lrate * ones(Dv));
        elseif strcmp(ratemethod,'componentwise_goal')                  
            actrate{j} = min(goaldiff / gradCompAvgs(1,j),lrate * ones(Dv));
        elseif strcmp(ratemethod,'only_goal')
            actrate{j} = goaldiff / gradCompAvgs(1,j);
        elseif strcmp(ratemethod,'only_rate')
            actrate{j} = lrate * ones(Dv);
        else
            exit('Invalid learning rate method: %s',ratemethod);
        end
        avgrates(1,j) = sum(sum(actrate{j}))/(Dv^2);                
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
    if ge.k>1
        subplot(ge.k,hor,hor+2);
        vsamp = reshape(samples(n,:,ge.k+1:ge.k+(ge.Dv*ge.B)),size(samples,2)*ge.B,ge.Dv);
        viewImage(cov(vsamp),'magnif',false);
        title('sample cov');
    end
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

function [grand_avg,comp_avgs] = componentAbsMeans(cc)
    k = size(cc,2);
    comp_avgs = zeros(1,k);
    for j=1:k
        comp_avgs(1,j) = comp_avgs(1,j) + mean(mean(abs(cc{j}),2),1);
    end
    grand_avg = mean(comp_avgs,2);  
end

function delPrint(num)
    for b=1:9+2*(floor(log10(num))+1)
        fprintf('\b');
    end
end

function ge = sortByGestalt(ge)
    indices = maxByRow(ge.G);
    newX = zeros(size(ge.X));
    %newV = zeros(size(ge.V));
    newG = zeros(size(ge.G));
    sofar = 0;
    for j = 1:ge.k
        gnum = sum(indices == j);
        newX(sofar + 1 : sofar + gnum,:,:) = ge.X(indices == j,:,:);
        %newv(sofar + 1 : sofar + gnum,:,:) = ge.V(indices == j,:,:);
        newG(sofar + 1 : sofar + gnum,:) = ge.G(indices == j,:);
        sofar = sofar + gnum;
    end
    ge.X = newX;
    %ge.V = newV;
    ge.G = newG;
end

function indices = maxByRow(A)
    indices = zeros(size(A,1),1);
    for i = 1:size(A,1)
        [~,indices(i,1)] = max(A(i,:));
    end
end