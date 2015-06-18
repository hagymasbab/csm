function [p,pexp] = stableMvnpdf(x,mu,C,varargin)
    parser = inputParser;   
    addParameter(parser,'logdetScaling','unknown');
    addParameter(parser,'scientific',false,@islogical);
    addParameter(parser,'invertedC',false,@islogical);
    addParameter(parser,'precompLogdet',NaN,@isnumeric);
    parse(parser,varargin{:}); 
    params = parser.Results; 

    if params.scientific && nargout < 2
        error('call with 2 output arguments to get scientific notation result');
    end
    d = size(x,1);
    
    if isnan(params.precompLogdet)
        ld = stableLogdet(C,'scaling',params.logdetScaling);
        if params.invertedC
            ld = -ld;
        end        
    else
        ld = params.precompLogdet;
    end    
    ld = -ld/2;
    
    v = x-mu;
    if params.invertedC
        quad = v' * C * v;
    elseif rcond(C) > 1e-16
        quad = v' * (C \ v);
    else
        quad = v' * pinv(C) * v;    
    end
    quad = -quad/2;    
    
    if params.scientific
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