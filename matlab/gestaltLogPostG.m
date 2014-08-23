function lp = gestaltLogPostG(g,V,ge,prior,precision)
    if (sum(g<0) > 0) || (strcmp(prior,'dirichlet') && ( (sum(g>1) > 0) || (sum(g) > 1) ))
        lp=-Inf;
        return
    end
    if size(g,1) == ge.k-1
        g = [g;1-sum(g)];
    end
    B = size(V,1);
    if ~precision
        Cv = componentSum(g,ge.cc);
        [~,err] = chol(Cv);
        if err > 0
            lp = -Inf;
            return;
        end
    else
        P = componentSum(g,ge.pc);
    end
    
    if strcmp(prior,'dirichlet')
        logprior = (ge.sparsity - 1) * sum(log(g));
    elseif strcmp(prior,'gamma')
        shape = 1;
        scale = 0.003;
        null_shape = 2;
        null_scale = 1;
        logprior = 0;
        for d = 1:ge.k            
            if ge.nullComponent && d == ge.k
                logprior = logprior + log(gampdf(g(d,1),null_shape,null_scale));
            else
                logprior = logprior + log(gampdf(g(d,1),shape,scale));
            end
        end
    end
    
    
    quad = 0;
    for b=1:B
        vb = V(b,:)';
        
        if ~precision
            quad = quad + vb' * (Cv \ vb);
            %this should be faster
%             opts.LT = true;
%             opts.UT = false;
%             temp = linsolve(U',vb,opts);
%             opts.LT = false;
%             opts.UT = true;
%             rightvec = linsolve(U,temp,opts);
%             quad = quad + vb' * rightvec;
        else
            quad = quad + vb' * P * vb;
        end
    end
    if ~precision
        lp = (-1/2) * ( B* log(det(Cv)) + quad ) + logprior;    
    else
        lp = (-1/2) * ( B* log(1/det(P)) + quad ) + logprior;
    end
end