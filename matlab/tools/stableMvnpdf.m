function p = stableMvnpdf(x,mu,C)
    d = size(x,1);
    ld = stableLogdet(C);
    v = x-mu;
    if log10(rcond(C)) < 10
        quad = v' * pinv(C) * v;
    else
        quad = v' * (C \ v);
    end
    exped = exp(-(ld + quad)/2);
    normconst = ((2*pi)^(-d/2));
    p = normconst * exped;
end