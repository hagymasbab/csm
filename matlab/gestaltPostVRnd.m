function V = gestaltPostVRnd(ge,xind,g)
    % construct the covariance and mean of the conditional posterior over v
    sAA = (1/ge.obsVar) * ge.AA;
    iCv = inv(componentSum(g,ge.cc));
    cov = inv(sAA + iCv);
    largecov = [];
    largemean = zeros(ge.B*ge.Dv,1);
    for b=1:ge.B
        largecov = blkdiag(largecov,cov);
        Ax = reshape(ge.tX(xind,b,:),1,ge.Dv)';
        % TEST: kivettem a minuszt
        m = (2/ge.obsVar) * cov * Ax;
        largemean(1+(b-1)*ge.Dv:b*ge.Dv,1) = m;
    end
    
    % generate a sample from this distribution
    V = mvnrnd(largemean',largecov);
    V = reshape(V,ge.B,ge.Dv);
end