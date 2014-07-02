function disc = testGestaltDer(ge)
    
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
            quad = quad + v(:,b)' * iC * v(:,b);
        end
        lp = -(b*log(det(C)) + quad) / 2;
    end

    function grad = loggaussgrad(C,v)
        b = size(v,2);
        iC = inv(C);
        quad = 0;
        for i=1:b
            quad = quad + iC * v(:,b) * v(:,b)' * iC;
        end
        grad = -(b*iC - quad) / 2;
    end

    function cv = kovmat(cholmat,g,ge)
        cv = cholmat2cov(cholmat,g,ge);
        cv = trace(cv);
    end

    function gradmat = kovmatder(cholmat,g,ge)
        gradcell = cell(1,ge.k);
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));
        for i=1:ge.k
            gradcell{i} = 2 * g(i,1) * choles{i};
        end
        gradmat = cell2mat(gradcell);
    end

    function lp = loggausschol(cholmat,g,v,ge)
        cv = cholmat2cov(cholmat,g,ge);
        lp = loggauss(cv,v);
    end

    function grad = chained(cholmat,g,v,ge)
        cv = cholmat2cov(cholmat,g,ge);
        cvgrad = loggaussgrad(cv,v);
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));
        gradcell = cell(1,ge.k);
        for i=1:ge.k
            gradcell{i} = 2 * g(i,1) * cvgrad * choles{i};
        end
        grad = cell2mat(gradcell);
    end

    function lp = gestaltUCDLL(cholmat,ge,samples)
        % unnormalised complete-data log-likelihood
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));        
        lp = gestaltCompleteDataLogLikelihood(ge,samples,choles);
    end

    function gradmat = gestaltDerUCDLL(cholmat,ge,samples)
        % derivative of unnormalised complete-data log-likelihood
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));
        grad = gestaltParamGrad(ge,samples,choles);
        gradmat = cell2mat(grad);
    end
    
    ge.B = 10;
    nSamp = 1;
    samples = zeros(1,nSamp,ge.k+ge.B*ge.Dv);
    samples(1,:,:) = gestaltGibbs(ge,1,nSamp);

    cholesky = cell(1,ge.k);
    for j=1:ge.k
        cholesky{j} = chol(ge.cc{j});
    end
    cholmat = cell2mat(cholesky);

    Cv = componentSum(0.5*ones(ge.k,1),ge.cc);
    v = mvnrnd(zeros(100,ge.Dv),Cv)';
    g = 0.5*ones(ge.k,1);
    
    a = @(x) gestaltUCDLL(x,ge,samples);
    b = @(x) gestaltDerUCDLL(x,ge,samples);
    init = cholmat;
%     a = @(x) kovmat(x,g,ge);
%     b = @(x) kovmatder(x,g,ge);

    a = @(x) loggausschol(x,g,v,ge);
    b = @(x) chained(x,g,v,ge);

%     a = @(x) loggauss(x,v);
%     b = @(x) loggaussgrad(x,v);
%     init = Cv;
    
    disc = checkDerivative(a,b,init,false);

end