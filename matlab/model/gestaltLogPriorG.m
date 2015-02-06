function lp = gestaltLogPriorG(g,ge)
    distribution = ge.prior;
    if strcmp(distribution,'gamma')
        lp = 0;
        for d = 1:ge.k            
            if ge.nullComponent && (d == ge.k)
                %lp = lp + log( gampdf(g(d,1),ge.null_shape,ge.null_scale) );
                lp = lp + (ge.null_shape-1) * log(g(d,1)) - g(d,1) / ge.null_scale;
            else
                %lp = lp + log( gampdf(g(d,1),ge.g_shape,ge.g_scale) );
                lp = lp + (ge.g_shape-1) * log(g(d,1)) - g(d,1) / ge.g_scale;
            end
        end
    elseif strcmp(distribution,'dirichlet')
         lp = (ge.sparsity - 1) * sum(log(g));
    else
        error('Prior distribution %s not implemented',distribution);
    end
end