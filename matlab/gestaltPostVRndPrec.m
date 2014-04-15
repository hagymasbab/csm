function V = gestaltPostVRndPrec(ge,xind,g)
    % construct the covariance and mean of the conditional posterior over v
    sAA = (1/ge.obsVar) * ge.AA;
    P = componentSum(g,ge.pc);
    cov = inv(sAA + P);
    %det(cov)
    
    V = zeros(ge.B,ge.Dv);
    for b=1:ge.B
        Ax = reshape(ge.tX(xind,b,:),1,ge.Dv)';
        m = (1/ge.obsVar) * cov * Ax;
        V(b,:) = mvnrnd(m',cov);
    end
end