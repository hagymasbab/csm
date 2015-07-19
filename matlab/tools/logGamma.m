function p = logGamma(x,alpha,theta)
    % up to a constant
    p = (alpha - 1) * sum(log(x)) - sum(x) / theta;
end