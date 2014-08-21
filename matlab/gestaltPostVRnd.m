function [V,mu] = gestaltPostVRnd(ge,xind,g,z,precision)

    % construct the covariance and mean of the conditional posterior over v
    sAA = ((z*z)/ge.obsVar) * ge.AA;
    if ~precision
        Cv = componentSum(g,ge.cc);
        cov = inv(sAA + inv(Cv));       
        % viewImage(cov);
    else
        P = componentSum(g,ge.pc);
        cov = inv(sAA + P);
    end
        
    V = zeros(ge.B,ge.Dv);
    for b=1:ge.B        
        ATx = ge.A' * reshape(ge.X(xind,b,:),ge.Dv,1);
        m = (z/ge.obsVar) * cov * ATx;
        V(b,:) = mvnrnd(m',cov);
    end
end