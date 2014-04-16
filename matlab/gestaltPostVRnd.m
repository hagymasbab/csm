function V = gestaltPostVRnd(ge,xind,g,precision)
    % construct the covariance and mean of the conditional posterior over v
    sAA = (1/ge.obsVar) * ge.AA;
    if ~precision
        P = inv(componentSum(g,ge.cc));
    else
        P = componentSum(g,ge.pc);
    end
    cov = inv(sAA + P);
    
    V = zeros(ge.B,ge.Dv);
    for b=1:ge.B
        Ax = reshape(ge.tX(xind,b,:),1,ge.Dv)';
        m = (1/ge.obsVar) * cov * Ax;
        V(b,:) = mvnrnd(m',cov);
    end
end