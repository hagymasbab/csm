function ld = stableLogdet(A,varargin)  
    parser = inputParser;   
    addParameter(parser,'scaling','unknown');
    parse(parser,varargin{:});       
    
    scale = 0;
    if strcmp(parser.Results.scaling,'up')
        scale = 1;
    elseif strcmp(parser.Results.scaling,'down')
        scale = -1;
    end
        
    if scale == 0
        dC = det(A); 
    end
    if scale == 0 && dC ~= 0 && abs(dC) ~= Inf
        ld = log(dC);
    else
        % if the determinant sucks, scale the fuck out of it
        orders = 700;
        d = size(A,1);        
        if scale == 1 || (scale == 0 && dC == 0) 
            s = 10^(orders/d);            
        elseif scale == -1 || (scale == 0 && abs(dC) == Inf)
            s = 10^(-orders/d);            
        end
        scaledDet = det(s*A);
        % if it still sucks, eigenvalues to the rescue
        if scaledDet == 0 || abs(scaledDet) == Inf                        
            ld = sum(log(eig(A)));
        else
            ld = log(scaledDet) - d * log(s); 
        end
    end
    
end
    