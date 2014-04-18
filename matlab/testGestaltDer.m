function disc = testGestaltDer(ge)
    
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

    samples = zeros(1,5,ge.k+ge.B*ge.Dv);
    samples(1,:,:) = gestaltGibbs(ge,1,5);

    cholesky = cell(1,ge.k);
    for j=1:ge.k
        cholesky{j} = chol(ge.cc{j});
    end
    cholmat = cell2mat(cholesky);

    a = @(x) gestaltUCDLL(x,ge,samples);
    b = @(x) gestaltDerUCDLL(x,ge,samples);
    
    disc = checkDerivative(a,b,cholmat);

end