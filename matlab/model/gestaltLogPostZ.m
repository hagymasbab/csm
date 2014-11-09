function lp = gestaltLogPostZ(z,xind,V,ge)
    if z < 0 
        lp = -Inf;
        return
    end
    
    if ndims(V) == 3
        V = reshape(V,ge.B,ge.Dv);
    end
    
    logprior = log( gampdf(z,ge.z_shape,ge.z_scale) );
    
    X = reshape(ge.X(xind,:,:),ge.B,ge.Dx);
    quad = 0;
    for b=1:ge.B
        vec = X(b,:)' - z * ge.A * V(b,:)';
        quad = quad + vec' * vec;
    end
    
    lp = (-1/2) * ( ge.B * ge.Dx * log(ge.obsVar) + quad / ge.obsVar) + logprior;
end