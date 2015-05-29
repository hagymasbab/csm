function [coefficient,exponent] = sciNot(a,logarithm)    
    s = sign(a);
    if logarithm
       % the original number cannot be negative, as we received the
       % logarithm of it
       a = a / log(10);       
       exponent = floor(a);
       coefficient = 10^(a - exponent);
    else        
        a = abs(a);
        a = log10(a);
        exponent = floor(a);
        coefficient = s .* 10.^(a-exponent);
    end
    
end