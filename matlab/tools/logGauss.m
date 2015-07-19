function lp = logGauss(x,mu,C)
    ld = stableLogdet(C);
    iC = stableInverse(C);
    lp = - (ld + (x - mu)' * iC * (x - mu))/2;
end