function ld = stableLogdet(A)
    
    dC = det(A);     
    if dC ~= 0 && abs(dC) ~= Inf
        ld = log(dC);
    else
        % if the determinant sucks, scale the fuck out of it
        orders = 700;
        d = size(A,1);        
        if dC == 0 
            s = 10^(orders/d);
        elseif abs(dC) == Inf
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
    