function lp = gestaltLogPostG(g,V,ge,precision)
    if (sum(g<0) > 0) || (sum(g>1) > 0) 
        lp=-Inf;
        return
    end
    if size(g,1) == ge.k-1
        g = [g;1-sum(g)];
    end
    B = size(V,1);
    if ~precision
        Cv = componentSum(g,ge.cc);
        [U,err] = chol(Cv);
        if err > 0
            lp = -Inf;
            return;
        end
    else
        P = componentSum(g,ge.pc);
    end
    prior = (ge.sparsity - 1) * sum(log(g));
    quad = 0;
    for b=1:B
        vb = V(b,:)';
        
        if ~precision
            %this should be faster
            opts.LT = true;
            opts.UT = false;
            temp = linsolve(U',vb,opts);
            opts.LT = false;
            opts.UT = true;
            rightvec = linsolve(U,temp,opts);
            quad = quad + vb' * rightvec;
        else
            quad = quad + vb' * P * vb;
        end
    end
    if ~precision
        lp = (-1/2) * ( B* log(det(Cv)) + quad ) + prior;    
    else
        lp = (-1/2) * ( B* log(1/det(P)) + quad ) + prior;
    end
end