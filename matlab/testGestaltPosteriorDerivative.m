function disc = testGestaltPosteriorDerivative(which,B)
    
    function lp = logpost_V(inputvec,ge,X)
        g = ge.G(1,:)';
        V = reshape(inputvec,ge.B,ge.Dv);
        z = ge.Z(1,1);
        iC = inv(componentSum(g,ge.cc));
        lp = gestaltFullLogPosterior(ge,X,V,g,z,iC);
    end

    function grad = logpostgrad_V(inputvec,ge,X)
        g = ge.G(1,:)';
        V = reshape(inputvec,ge.B,ge.Dv);
        z = ge.Z(1,1);
        iC = inv(componentSum(g,ge.cc));
        grad = gestaltFullLogPosteriorGrad(ge,X,V,g,z,iC);
        grad = grad(ge.k+1:end-1,1);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function lp = logpost_G(inputvec,ge,X)
        g = inputvec;
        V = reshape(ge.V(1,:,:),ge.B,ge.Dv);
        z = ge.Z(1,1);
        iC = inv(componentSum(g,ge.cc));
        lp = gestaltFullLogPosterior(ge,X,V,g,z,iC);
    end

    function grad = logpostgrad_G(inputvec,ge,X)
        g = inputvec;
        V = reshape(ge.V(1,:,:),ge.B,ge.Dv);
        z = ge.Z(1,1);
        iC = inv(componentSum(g,ge.cc));
        grad = gestaltFullLogPosteriorGrad(ge,X,V,g,z,iC);
        grad = grad(1:ge.k,1);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function lp = logpost_Z(inputvec,ge,X)
        g = ge.G(1,:)';
        V = reshape(ge.V(1,:,:),ge.B,ge.Dv);
        z = inputvec;
        iC = inv(componentSum(g,ge.cc));
        lp = gestaltFullLogPosterior(ge,X,V,g,z,iC);
    end

    function grad = logpostgrad_Z(inputvec,ge,X)
        g = ge.G(1,:)';
        V = reshape(ge.V(1,:,:),ge.B,ge.Dv);
        z = inputvec;
        iC = inv(componentSum(g,ge.cc));
        grad = gestaltFullLogPosteriorGrad(ge,X,V,g,z,iC);
        grad = grad(end,1);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function lp = logpost(inputvec,ge,X)
        g = inputvec(1:ge.k,1);
        V = reshape(inputvec(ge.k+1:end-1,1),ge.B,ge.Dv);
        z = inputvec(end,1);
        iC = inv(componentSum(g,ge.cc));
        lp = gestaltFullLogPosterior(ge,X,V,g,z,iC);
    end

    function grad = logpostgrad(inputvec,ge,X)
        g = inputvec(1:ge.k,1);
        V = reshape(inputvec(ge.k+1:end-1,1),ge.B,ge.Dv);
        z = inputvec(end,1);
        iC = inv(componentSum(g,ge.cc));
        grad = gestaltFullLogPosteriorGrad(ge,X,V,g,z,iC);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define model
    ge = gestaltCreate('temp','Dx',256,'k',3,'nullComponent',true,'generate',true,'B',B,'N',1);    
    
    % define stimulus
    X = reshape(ge.X(1,:,:),ge.B,ge.Dx);
    
    if strcmp(which,'full')    
        inputvec_size = ge.k+ge.B*ge.Dv+1;        
        a = @(inputvec) logpost(inputvec,ge,X);
        b = @(inputvec) logpostgrad(inputvec,ge,X);
    elseif strcmp(which,'g')    
        inputvec_size = ge.k;
        a = @(inputvec) logpost_G(inputvec,ge,X);
        b = @(inputvec) logpostgrad_G(inputvec,ge,X);        
    elseif strcmp(which,'v')    
        inputvec_size = ge.B*ge.Dv;
        a = @(inputvec) logpost_V(inputvec,ge,X);
        b = @(inputvec) logpostgrad_V(inputvec,ge,X);
    elseif strcmp(which,'z')    
        inputvec_size = 1;
        a = @(inputvec) logpost_Z(inputvec,ge,X);
        b = @(inputvec) logpostgrad_Z(inputvec,ge,X);
    end
    
    init = ones(inputvec_size,1) * 0.1;        
    disc = checkDerivative(a,b,init,false);
    
end