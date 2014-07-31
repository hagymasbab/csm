function [V,mu] = gestaltPostVRnd(ge,xind,g,precision)
    % there could be an approximation, but it's assumptions are not met, so
    % it doesn't work
    approx = false;
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
        %Ax = reshape(ge.tX(xind,b,:),1,ge.Dv)';
        ATx = ge.A' * reshape(ge.X(xind,b,:),ge.Dv,1);
        m = (1/ge.obsVar) * cov * ATx;
        V(b,:) = mvnrnd(m',cov);
    end
end