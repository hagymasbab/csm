function R = gaussianRandom(N,m,precision)
    means = repmat(m',N,1);
    dim = size(m,1);
    standards = normrnd(0,1,[dim,N]);
    L = chol(precision,'lower');
    R = means + (L * standards)';
end