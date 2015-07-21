function grad = gestaltLogPostGZGrad(g,z,x,ge)
    
    Cv = componentSum(g,ge.cc);
    ACAT = ge.A * Cv * ge.A';
    Cx = ge.obsVar * eye(ge.Dv) + z^2 * ACAT;    
    
    iCx = stableInverse(Cx);
    leftmat = iCx - iCx * (x * x') * iCx;

    % derivative w.r.t. g    
    ggrad = zeros(ge.k,1);
    for kk=1:ge.k
        right_gk = ge.A * ge.cc{kk} * ge.A';
        tr = trace(leftmat*right_gk);
        ggrad(kk) = - (z^2 / 2) * tr;
    end
    ggrad = ggrad + logGammaDerivative(g,ge.g_shape,ge.g_scale);
    
    % derivative w.r.t. z
    tr = trace(leftmat*ACAT);
    zgrad = -z * tr + logGammaDerivative(z,ge.z_shape,ge.z_scale);
    
    grad = [ggrad;zgrad];
end
        