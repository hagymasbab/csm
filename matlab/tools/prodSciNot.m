function [coefficient,exponent] = prodSciNot(a,e)
    exponent = 0;
    coefficient = 1;
    for i=1:size(a,2)
        coefficient = coefficient * a(i);
        if coefficient >= 10
            exponent = exponent + 1;
            coefficient = coefficient / 10;
        elseif coefficient <= -10
            exponent = exponent - 1;
            coefficient = coefficient / 10;
        end
    end
    exponent = exponent + sum(e,2);
end