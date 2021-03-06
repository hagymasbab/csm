function disc = testGestaltDer(ge,formula,randseed)
    
    function cv = cholmat2cov(cholmat,g,ge)
        cv = zeros(ge.Dv);
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));
        for i=1:ge.k
            cv = cv + g(i,1) * (choles{i}' * choles{i});
        end
    end

    function lp = loggauss(C,v)
        b = size(v,2);
        iC = inv(C);
        quad = 0;
        for i=1:b
            v_act = v(:,i);
            %size(v_act)
            quad = quad + v_act' * iC * v_act;
        end
        lp = -(b*log(det(C)) + quad) / 2;
    end

    function grad = loggaussgrad(C,v)
        b = size(v,2);
        iC = inv(C);
        quad = zeros(size(C));
        for i=1:b
            v_act = v(:,i);
            quad = quad + iC * (v_act * v_act') * iC;
        end
        grad = -(b*iC - quad) / 2;
    end

    function cv = kovmat(upperleft,cholmat,g,ge)        
        cholmat(1,1) = upperleft;
        cv = cholmat2cov(cholmat,g,ge);        
        %cv = trace(cv);
    end

    function gradmat = kovmatder(upperleft,cholmat,g,ge)        
        gradcell = cell(1,ge.k);
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));
        choles{1}(1,1) = upperleft;
        U_hat = derivQuadByElement(choles{1},1,1);
        gradmat = g(1) * U_hat;
%         for i=1:ge.k
%             gradcell{i} = 2 * g(i,1) * choles{i};
%         end
%         gradmat = cell2mat(gradcell);
    end

    function lp = loggausschol(upperleft,cholmat,g,v,ge)
        cholmat(1,1) = upperleft;
        cv = cholmat2cov(cholmat,g,ge);
        %cv = cholmat;
        lp = loggauss(cv,v);
    end

    function grad = chained(upperleft,cholmat,g,v,ge)
        % substituting th input value
        cholmat(1,1) = upperleft;
        % assembling the covariance matrix
        cv = cholmat2cov(cholmat,g,ge);
        %gradient of the log-gaussian formula w.r.t the covariance matrix
        dLdC = loggaussgrad(cv,v);
                
        % cell representation of the covariance components
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));
        % derivative of U^TU w.r.t. the first element of U_1
        U_hat = derivQuadByElement(choles{1},1,1);
        % derivative of the covariance matrix w.r.t. the input
        dCdu = g(1,1) * U_hat;
        
        % derivative of the log-gaussian w.r.t. the input
        grad = trace(dLdC' * dCdu);
    end

    function lp = gestaltUCDLL(upperleft,cholmat,ge,samples,precision)
        % unnormalised complete-data log-likelihood
        cholmat(1,1) = upperleft;
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));        
        lp = gestaltCompleteDataLogLikelihood(ge,samples,choles,'precision',precision);
    end

    function elementGrad = gestaltDerUCDLL(upperleft,cholmat,ge,vsamp,gsamp,precision)
        % derivative of unnormalised complete-data log-likelihood
        cholmat(1,1) = upperleft;
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));
        grad = gestaltParamGrad(ge,vsamp,gsamp,choles,'precision',precision);        
        elementGrad = grad{1}(1,1);
        %gradmat = cell2mat(grad);
    end

    if strcmp(randseed,'last')
        load lastrandseed;
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');
    
    nSamp = 10;
    samples = zeros(1,nSamp,ge.k+ge.B*ge.Dv);
    [vsamp,gsamp,zsamp] = gestaltGibbs(ge,1,nSamp);
    s = mergeSamples(vsamp,gsamp,zsamp);
    samples(1,:,:) = s(:,1:end-1);

    cholesky = cell(1,ge.k);
    for ak=1:ge.k
        cholesky{ak} = chol(ge.cc{ak});
    end
    cholmat = cell2mat(cholesky);

    Cv = componentSum(0.5*ones(ge.k,1),ge.cc);
    v = mvnrnd(zeros(ge.B,ge.Dv),Cv)';
    g = 0.5*ones(ge.k,1);
    %samples(1,1,:) = [g' reshape(v,1,ge.B*ge.Dv)];

    init = 1;
    
    if strcmp(formula,'cdll')    
        % R -> R
        a = @(x) gestaltUCDLL(x,cholmat,ge,samples,false);
        b = @(x) gestaltDerUCDLL(x,cholmat,ge,vsamp,gsamp,false);        
        
    elseif strcmp(formula,'cdll_prec')    
        % R -> R
        a = @(x) gestaltUCDLL(x,cholmat,ge,samples,true);
        b = @(x) gestaltDerUCDLL(x,cholmat,ge,samples,true);        
        
    elseif strcmp(formula,'kovmat')
        % R -> R^N
        a = @(x) kovmat(x,cholmat,g,ge);
        b = @(x) kovmatder(x,cholmat,g,ge);        
        
    elseif strcmp(formula,'chain')
        % R -> R
        a = @(x) loggausschol(x,cholmat,g,v,ge);
        b = @(x) chained(x,cholmat,g,v,ge);        
        
    elseif strcmp(formula,'loggauss')
        % R^N -> R
        a = @(x) loggauss(x,v);
        b = @(x) loggaussgrad(x,v);
        init = Cv;
    else
        error('not valid formula');
    end
    
    disc = checkDerivative(a,b,init,false);

end