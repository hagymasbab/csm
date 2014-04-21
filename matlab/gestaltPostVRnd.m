function V = gestaltPostVRnd(ge,xind,g,precision,approx)
    % construct the covariance and mean of the conditional posterior over v
    sAA = (1/ge.obsVar) * ge.AA;
    if ~precision
        Cv = componentSum(g,ge.cc);
        if ~approx
            cov = inv(sAA + inv(Cv));
        else
            cov = Cv - Cv * sAA * Cv;
        end
    else
        P = componentSum(g,ge.pc);
        cov = inv(sAA + P);
    end
    
    
    V = zeros(ge.B,ge.Dv);
    for b=1:ge.B
        Ax = reshape(ge.tX(xind,b,:),1,ge.Dv)';
        m = (1/ge.obsVar) * cov * Ax;
        V(b,:) = mvnrnd(m',cov);
    end
end