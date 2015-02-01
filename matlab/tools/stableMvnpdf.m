function [p,pexp] = stableMvnpdf(x,mu,C,scientific)
    if scientific && nargout < 2
        error('call with 2 output arguments to get scientific notation result');
    end
    d = size(x,1);
    ld = stableLogdet(C);
    ld = -ld/2;
    v = x-mu;
    rc = rcond(C);
    if log10(rc) < -30
        quad = v' * pinv(C) * v;
    else
        quad = v' * (C \ v);
    end
    quad = -quad/2;
    normconst = ((2*pi)^(-d/2));
    
    if scientific
        [ld_coeff,ld_exp] = sciNot(ld,true);
        [quad_coeff,quad_exp] = sciNot(quad,true);
        [norm_coeff,norm_exp] = sciNot(normconst,false);
        [p,pexp] = prodSciNot([ld_coeff quad_coeff norm_coeff],[ld_exp quad_exp norm_exp]);
    else
        exped = exp(ld + quad);          
        p = normconst * exped;  
    end
      
    
end