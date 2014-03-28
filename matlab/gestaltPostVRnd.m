function v = gestaltPostVRnd(ge,xind,g)
    % construct the covariance and mean of the conditional posterior over v
    sAA = (1/ge.obsVar) * ge.AA;
    Ax = ge.tX(xind,:)';
    iCv = inv(componentSum(g,ge.cc));
    cov = inv(sAA + iCv);
    % TEST: kivettem a minuszt
    m = (2/ge.obsVar) * cov * Ax;
    % generate a sample from this distribution
    v = mvnrnd(m',cov)';
end