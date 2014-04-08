function lp = gestaltLogPostG(g,V,ge)
    if (sum(g<0) > 0) || (sum(g>1) > 0) 
        lp=-Inf;
        return
    end
    if size(g,1) == ge.k-1
        g = [g;1-sum(g)];
    end
    B = size(V,1);
    Cv = componentSum(g,ge.cc);
    %iCv = inv(Cv);
    iCv = eye(ge.Dv) / Cv;
    prior = (ge.sparsity - 1) * sum(log(g));
    quad = 0;
    for b=1:B
        vb = V(b,:)';
        %quad = quad + vb'*(Cv \ vb);
        quad = quad + vb'* iCv * vb;
    end
    lp = (-1/2) * ( B* log(abs(det(Cv))) + quad ) + prior;
    %lp = (-1/2) * ( B* log(abs(det(Cv))) + quad );
    %lp = prior;
end