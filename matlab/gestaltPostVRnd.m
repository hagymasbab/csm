function V = gestaltPostVRnd(ge,xind,g)
    % construct the covariance and mean of the conditional posterior over v
    sAA = (1/ge.obsVar) * ge.AA;
    iCv = inv(componentSum(g,ge.cc));
    prec = sAA + iCv;
    cov = inv(sAA + iCv);
    %cholprec = chol(prec,'lower');    
    
    V = zeros(ge.B,ge.Dv);
    for b=1:ge.B
        Ax = reshape(ge.tX(xind,b,:),1,ge.Dv)';
        m = (1/ge.obsVar) * cov * Ax;
        V(b,:) = mvnrnd(m',cov);
    end
end