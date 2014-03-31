function lp = gestaltLogPostG(g,v,ge)
    if (sum(g<0) > 0) || (sum(g>1) > 0) 
        lp=-Inf;
        return
    end
    if size(g,1) == ge.k-1
        g = [g;1-sum(g)];
    end
    Cv = componentSum(g,ge.cc);
    prior = (ge.sparsity - 1) * sum(log(g));
    lp = (-1/2) * ( log(abs(det(Cv))) + v'*(Cv \ v) ) + prior;
    %lp = (-1/2) * ( log(abs(det(Cv))) + v'*(Cv \ v) );
end