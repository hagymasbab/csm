function V = gestaltPostVRnd(ge,xind,g)
    % construct the covariance and mean of the conditional posterior over v
    sAA = (1/ge.obsVar) * ge.AA;
    iCv = inv(componentSum(g,ge.cc));
    cov = inv(sAA + iCv);
    
    V = zeros(ge.B,ge.Dv);
    for b=1:ge.B
        Ax = reshape(ge.tX(xind,b,:),1,ge.Dv)';
        % TEST: kivettem a minuszt
        m = (2/ge.obsVar) * cov * Ax;
        V(b,:) = mvnrnd(m',cov);
    end
end