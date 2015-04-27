function [V,m] = gestaltPostVRnd(ge,xind,g_or_Cv,z,precision)

    % construct the covariance and mean of the conditional posterior over v
    sAA = ((z*z)/ge.obsVar) * ge.AA;
    if ~precision
        if size(g_or_Cv,2) == ge.Dv
            Cv = g_or_Cv;
        else
            Cv = componentSum(g_or_Cv,ge.cc);        
        end
        %icv = Cv \ eye(ge.Dv);
        %icovm = sAA + icv;
        %icovm = sAA + inv(Cv);
        %covm = icovm \ eye(ge.Dv);
        icv = inv(Cv);
        covm = inv(sAA + icv);       
        % viewImage(cov);
    else
        % TODO implement the g_or_Cv thing
        P = componentSum(g,ge.pc);
        covm = inv(sAA + P);
    end
        
    V = zeros(ge.B,ge.Dv);
    for b=1:ge.B        
        ATx = ge.A' * reshape(ge.X(xind,b,:),ge.Dv,1);
        m = (z/ge.obsVar) * covm * ATx;
        V(b,:) = mvnrnd(m',covm);
    end
end