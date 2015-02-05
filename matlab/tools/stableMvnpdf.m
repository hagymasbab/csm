function [p,pexp] = stableMvnpdf(x,mu,C,scientific)
    if scientific && nargout < 2
        error('call with 2 output arguments to get scientific notation result');
    end
    d = size(x,1);
    ld = stableLogdet(C);
    ld = -ld/2;
    v = x-mu;
    rc = rcond(C);
    if log10(rc) < -15
        quad = v' * pinv(C) * v;
    else
        quad = v' * (C \ v);
    end
    quad = -quad/2;    
    
    if scientific
        normconst = (-d/2) * log(2*pi);
        [norm_coeff,norm_exp] = sciNot(normconst,true);
        
        [ld_coeff,ld_exp] = sciNot(ld,true);
        [quad_coeff,quad_exp] = sciNot(quad,true);        
        [p,pexp] = prodSciNot([ld_coeff quad_coeff norm_coeff],[ld_exp quad_exp norm_exp]);        
    else
        normconst = ((2*pi)^(-d/2));
        exped = exp(ld + quad);          
        p = normconst * exped;  
    end
      
    
end