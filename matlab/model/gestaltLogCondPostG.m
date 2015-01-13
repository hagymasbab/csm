function lcp = gestaltLogCondPostG(gi,g,compIdx,V,ge,prior,precision,cv,cc)
    if (gi < 0 || sum(g<0) > 0)
        lcp=-Inf;
        return
    end
    
    actdiff = g(compIdx,1) - gi;
    cv = cv - actdiff * cc{compIdx};
    g(compIdx,1) = gi;

    loglike = gestaltLogLikeV(V,g,ge,precision,[],cv);
    
    if strcmp(prior,'gamma')
        if ge.nullComponent && (compIdx == ge.k)
            lp = log( gampdf(gi,ge.null_shape,ge.null_scale) );
        else
            lp = log( gampdf(gi,ge.g_shape,ge.g_scale) );
        end
    else
        throw(MException('Gestalt:LogCondPost:NotImplemented','Log-conditional posteriror for g is only implemented for Gamma priors.'));
    end
    
    lcp = loglike + lp;
    %lcp = loglike;
    %lcp = lp;
end