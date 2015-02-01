function [coefficient,exponent] = sciNot(a,logarithm)    
    s = sign(a);
    if logarithm
       a = a / log(10);       
       exponent = floor(a);
       coefficient = 10^(a+s*exponent);
    else        
        a = abs(a);
        a = log10(a);
        exponent = floor(a);
        coefficient = s * 10^(a-exponent);
    end
    
end