function lp = gestaltLogPostGPrec(g,V,ge)
    if (sum(g<0) > 0) || (sum(g>1) > 0) 
        lp=-Inf;
        return
    end
    if size(g,1) == ge.k-1
        g = [g;1-sum(g)];
    end
    B = size(V,1);
    P = componentSum(g,ge.pc);
    
    prior = (ge.sparsity - 1) * sum(log(g));
    quad = 0;
    for b=1:B
        vb = V(b,:)';
        quad = quad + vb' * P * vb;
    end
    lp = (-1/2) * ( B* log(abs(1/det(P))) + quad ) + prior;
end