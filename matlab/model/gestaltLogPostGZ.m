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
    ld = stableLogdet(Cx);
    iCx = stableInverse(Cx);
%     [~,e] = cholcov(Cx)
%         [~,e] = cholcov(iCx)

    likelihood = - (ld + x' * iCx * x)/2;
    gprior = (ge.g_shape - 1) * sum(log(g)) - sum(g) / ge.g_scale;
    zprior = (ge.z_shape - 1) * log(z) - z / ge.z_scale;
    lp = likelihood + gprior + zprior;
end