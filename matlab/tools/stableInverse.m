function iA = stableInverse(A)
    if log10(rcond(A)) < -16
        iA = pinv(A);
    else
        % inv is by no means numerically stable
        %iA = inv(A);
        iA = eye(size(A,1)) / A;
        % we might mess up symmetry with mldivide for some reason, so
        if issymmetric(A)
            iA = 0.5 * (iA + iA');
        end
    end
end