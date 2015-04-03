function lp = gestaltFullLogPosterior(ge,X,V,g,z,iC)
    % no precision formulation is available
    % likelihood of X
    scalar = 0;
    for b = 1: ge.B
        zAv = z*ge.A*V(b,:)';
        scalar = scalar + (X(b,:)' - zAv)' * (X(b,:)' - zAv);
    end
    lx = - scalar / (2 * ge.obsVar);
    % likelihood of V
    lv = gestaltLogLikeV(V,g,ge,false,iC,[]);
    % prior of g
    lg = gestaltLogPriorG(g,ge);
    % prior of z
    lz = gestaltLogPriorZ(z,ge);
    
    lp = lx + lv + lg + lz;
end