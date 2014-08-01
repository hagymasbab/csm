function [V,mu] = gestaltPostVRnd(ge,xind,g,precision)

    % construct the covariance and mean of the conditional posterior over v
    sAA = (1/ge.obsVar) * ge.AA;
    if ~precision
        Cv = componentSum(g,ge.cc);
        cov = inv(sAA + inv(Cv));       
        viewImage(cov);
    else
        P = componentSum(g,ge.pc);
        cov = inv(sAA + P);
    end
        
    V = zeros(ge.B,ge.Dv);
    for b=1:ge.B        
        ATx = ge.A' * reshape(ge.X(xind,b,:),ge.Dv,1);
        m = (1/ge.obsVar) * cov * ATx;
        V(b,:) = mvnrnd(m',cov);
    end
end