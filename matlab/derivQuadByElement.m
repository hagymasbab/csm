function X_hat = derivQuadByElement(X,i,j)
    JJ = zeros(size(X));
    JJ(i,j) = 1;
    X_hat = X' * JJ + JJ' * X;
end