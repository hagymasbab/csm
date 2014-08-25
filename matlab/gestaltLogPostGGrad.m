function grad = gestaltLogPostGGrad(g,ge,V)
    Cv = componentSum(g,ge.cc);
    iCv = inv(Cv);
    quad = zeros(ge.Dv);
    for b=1:ge.B
        v = V(b,:)';
        quad = quad + iCv * (v * v') * iCv;
    end
    dldC = ge.B * iCv - quad;
    
    grad = zeros(size(g));
    for kk = 1:ge.k
        if ge.nullComponent && kk = ge.k
            shape = ge.null_shape;
            scale = ge.null_scale;
        else
            shape = ge.g_shape;
            scale = ge.g_scale;
        end
        act_g = g(kk,1);
        dpdg = (1/(gamma(shape)*scale^shape)) * exp(-act_g/scale) * ( act_g^(shape-2)*(shape-1) - act_g^(shape-1)/scale );
        grad(kk,1) = (-1/2) * trace(dldC*ge.cc{kk}) + dpdg;
    end
end