function iA = stableInverse(A)
    if log10(rcond(A)) < -15
        iA = pinv(A);
    else
        iA = inv(A);
    end
end