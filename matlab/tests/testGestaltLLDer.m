function testGestaltLLDer(formula,randseed,like_method)
    
    function lp = gestaltLLX(sigma_x,choles,ge,L,lm)
        ge.obsVar = sigma_x;
        lp = gestaltLogLikelihood2(ge,L,ge.X,choles,'loadSamples',true,'method',lm);
    end

    function grad = gestaltDerLLX(sigma_x,choles,ge,L)
        % derivative of log-likelihood
        grad = gestaltLogLikelihoodGradX(ge,L,ge.X,choles,sigma_x,'loadSamples',true);
    end

    function lp = gestaltLL(x,cholmat,ge,L,k,c1,c2,lm)
        % log-likelihood        
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));  
        choles{k}(c1,c2) = x;
        %lp = gestaltLogLikelihood(ge,L,ge.X,'cholesky',choles,'loadSamples',false);
        lp = gestaltLogLikelihood2(ge,L,ge.X,choles,'loadSamples',true,'method',lm);
    end

    function grad_el = gestaltDerLL(x,cholmat,ge,L,k,c1,c2)
        % derivative of log-likelihood
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));
        choles{k}(c1,c2) = x;
        grad = gestaltLogLikelihoodGradient(ge,L,ge.X,choles,'loadSamples',true);
        grad_el = grad{k}(c1,c2);
    end

    function p = gaussian(upperleft,kovmat,x,mu,c1,c2)
        kovmat(c1,c2) = upperleft;
        kovmat(c2,c1) = upperleft;
        p = mvnpdf(x',mu',kovmat);
    end

    function gradElement = gaussGrad(upperleft,kovmat,x,mu,c1,c2)
        kovmat(c1,c2) = upperleft;
        kovmat(c2,c1) = upperleft;
        kovmat
        p = mvnpdf(x',mu',kovmat);
        iC = inv(kovmat);     
        t = iC*(x-mu);
        grad = -p/2 * (iC - t*t');
        gradElement = grad(c1,c2);
    end

    setrandseed(randseed);
    ge = gestaltCreate('temp','Dx',64,'k',2,'N',10,'filters','OF','obsVar',0.5,'g_shape',1,'g_scale',1,'z_shape',2,'z_scale',2,'generateComponents',true,'generateData',true);
    U = cellchol(ge.cc);    
    L = 20;
    
    cholmat = cell2mat(U);
    
    
    if strcmp(formula,'ll')
        
        k = 1;
        c1 = 3;
        c2 = 4;
        % just to save the right samples
        lp = gestaltLogLikelihood2(ge,L,ge.X,U,'loadSamples',false);
        a = @(x) gestaltLL(x,cholmat,ge,L,k,c1,c2,like_method);
        b = @(x) gestaltDerLL(x,cholmat,ge,L,k,c1,c2);   
        
        init = 0.1;
        
    elseif strcmp(formula,'llx')
        
        lp = gestaltLogLikelihood2(ge,L,ge.X,U,'loadSamples',false);
        a = @(x) gestaltLLX(x,U,ge,L,like_method);
        b = @(x) gestaltDerLLX(x,U,ge,L);   
        
        init = 0.1;
        
    elseif strcmp(formula,'gausscov')
            
%         %kovmat = componentSum(1,ge.cc);
%         rm = randn(ge.Dv);
%         %kovmat = rm*rm';
%         kovmat = eye(ge.Dv) + 0.01 * (rm* rm');
%         y = squeeze(ge.X(1,1,:));
%         mu = zeros(ge.Dx,1);
        
        kovmat = [30 0.5;0.5 40];
        mu = [0;0];
        y = [1;1];
    
        a = @(x) gaussian(x,kovmat,y,mu,1,2);
        b = @(x) gaussGrad(x,kovmat,y,mu,1,2);   
    
        init = 0;
    
    else
        
        error('no such formula');
        
    end
    
    disc = checkDerivative(a,b,init,false)

end