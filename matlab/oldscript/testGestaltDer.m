function disc = testGestaltDer(cc)

    function lp = gestaltUCDLL(g,v,cc1,cc2)
        % unnormalised complete-data log-likelihood
        ccc{1} = cc1;
        ccc{2} = cc2;
        cv = componentSum(g,ccc);
        lp = (-1/2) * ( logdet(cv,'chol') + v'*inv(cv)*v);
    end

    function lp = gestaltDerUCDLL(g,v,cc1,cc2,j)
        % derivative of unnormalised complete-data log-likelihood
        ccc{1} = cc1;
        ccc{2} = cc2;
        icv = inv(componentSum(g,ccc));
        lp = (-1/2) * g(j) * (icv - icv * (v*v') * icv );
    end
    
    k = size(cc,2);
    d = size(cc{1},1);
    g = 0.5*ones(k,1);
    v = 0.5*ones(d,1);
    x = cc{1};
    cc2 = cc{2};

    a = @(x) gestaltUCDLL(g,v,x,cc2);
    b = @(x) gestaltDerUCDLL(g,v,x,cc2,1);
    
    disc = checkDerivative(a,b,x);

end