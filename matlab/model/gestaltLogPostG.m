function lp = gestaltLogPostG(g,V,ge,prior,precision)
    % log-posterior of g up to a constant
    
    if (sum(g<0) > 0) || (strcmp(prior,'dirichlet') && ( (sum(g>1) > 0) || (sum(g) > 1) ))
        lp=-Inf;
        return
    end
    
    if size(g,1) == ge.k-1 && strcmp(prior,'dirichlet')
        g = [g;1-sum(g)];
    end
    
    logprior = gestaltLogPriorG(g,ge,prior);        
    
    loglike = gestaltLogLikeV(V,g,ge,precision);
    
    lp = loglike + logprior;
end