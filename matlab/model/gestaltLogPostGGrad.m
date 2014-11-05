function grad = gestaltLogPostGGrad(g,V,ge,prior,precision)
    if ndims(V) == 3
        V = reshape(V,ge.B,ge.Dv);
    end
    if ~precision
        Cv = componentSum(g,ge.cc);
        iCv = inv(Cv);
    
        quad = zeros(ge.Dv);
        for b=1:ge.B
            v = V(b,:)';
            quad = quad + iCv * (v * v') * iCv;
        end
        dldC = ge.B * iCv - quad;
    else
        exit('Gradient of log-posterior over g not implemented for the precision formulation.');
    end
    
    grad = zeros(size(g));
    for kk = 1:ge.k
        act_g = g(kk,1);
        if strcmp(prior,'gamma')
            if ge.nullComponent && kk == ge.k
                shape = ge.null_shape;
                scale = ge.null_scale;
            else
                shape = ge.g_shape;
                scale = ge.g_scale;
            end        
            dpdg = (1/(gamma(shape)*scale^shape)) * exp(-act_g/scale) * ( act_g^(shape-2)*(shape-1) - act_g^(shape-1)/scale );
        else
            exit('Gradient of log-posterior over g not implemented for the specified prior.');
        end
        grad(kk,1) = (-1/2) * trace(dldC*ge.cc{kk}) + dpdg;
    end
end