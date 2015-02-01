function [coefficient,exponent] = sumSciNot(c1,e1,c2,e2)
    if c1 == 0 && e1 == 0
        coefficient = c2;
        exponent = e2;
        return;
    elseif c2 == 0 && e2 == 0
        coefficient = c1;
        exponent = e1;
        return;
    end
    exponent = max(e1,e2);
    expdiff = abs(e1-e2);
    if(e1 > e2)
        coefficient = c1 + 10^(-expdiff) * c2;
    else
        coefficient = c2 + 10^(-expdiff) * c1;
    end
    [coefficient,ex] = sciNot(coefficient,false);
    exponent = exponent + ex;
end
