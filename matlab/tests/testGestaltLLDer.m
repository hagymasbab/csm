function testGestaltLLDer(formula,randseed)
    
    function lp = gestaltLL(upperleft,cholmat,ge,L)
        % log-likelihood
        cholmat(1,1) = upperleft;
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));        
        lp = gestaltLogLikelihood(ge,L,ge.X,'cholesky',choles,'loadSamples',false);
    end

    function elementGrad = gestaltDerLL(upperleft,cholmat,ge,L)
        % derivative of log-likelihood
        cholmat(1,1) = upperleft;
        choles = mat2cell(cholmat,ge.Dv,ge.Dv*ones(1,ge.k));
        grad = gestaltLogLikelihoodGradient(ge,L,ge.X,choles,'loadSamples',true);
        elementGrad = grad{1}(1,1);        
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
    ge = gestaltCreate('temp','Dx',64,'k',10,'N',50,'filters','OF','obsVar',0.5,'g_shape',1,'g_scale',1,'z_shape',2,'z_scale',2,'generateComponents',true,'generateData',true);
    U = cellchol(ge.cc);    
    L = 10;
    
    cholmat = cell2mat(U);
    
    
    if strcmp(formula,'ll')

        a = @(x) gestaltLL(x,cholmat,ge,L);
        b = @(x) gestaltDerLL(x,cholmat,ge,L);   
        
        init = 1;
        
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
    end
    
    disc = checkDerivative(a,b,init,false)

end