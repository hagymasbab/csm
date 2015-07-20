function dervec = logGammaDerivative(x,alpha,theta)
    dervec = (alpha - 1) ./ x - 1 / theta;
end
