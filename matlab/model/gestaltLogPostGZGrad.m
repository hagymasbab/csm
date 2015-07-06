function grad = gestaltLogPostGZGrad(g,z,x,ge)
    
    Cv = componentSum(g,ge.cc);
    ACAT = ge.A * Cv * ge.A';
    Cx = ge.obsVar * eye(ge.Dv) + z^2 * ACAT;    
    
    iCx = stableInverse(Cx);
    leftmat = iCx - iCx * (x * x') * iCx;
%     leftmat = Cx \ (eye(ge.Dv) - (x * x') / Cx);

    % derivative w.r.t. g    
    ggrad = zeros(ge.k,1);
    for kk=1:ge.k
        right_gk = ge.A * ge.cc{kk} * ge.A';
        tr = trace(leftmat*right_gk);
        ggrad(kk) = - (z^2 / 2) * tr + (ge.g_shape-1)/g(kk) - 1/ge.g_scale;
    end
    
    % derivative w.r.t. z
    tr = trace(leftmat*ACAT);
    zgrad = -z * tr + (ge.z_shape-1)/z - 1/ge.z_scale;
    
    grad = [ggrad;zgrad];
end
        