function lp = gestaltLogPostGZ(g,z,x,ge)
    if ~strcmp(ge.prior,'gamma') || ge.B>1
        error('not implemented');
    end
    if any(g<0) || z<0
        lp = -Inf;
        return;
    end
    Cv = componentSum(g,ge.cc);
    Cx = ge.obsVar * eye(ge.Dx) + z^2 * ge.A * Cv * ge.A';
    lp = logGauss(x,0,Cx) + logGamma(g,ge.g_shape,ge.g_scale) + logGamma(z,ge.z_shape,ge.z_scale);    
end