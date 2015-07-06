function iA = stableInverse(A)
    if rcond(A) < 1e-16
        iA = pinv(A);
    else
        % inv is by no means numerically stable
        %iA = inv(A);
        iA = eye(size(A,1)) / A;        
    end
    % we might mess up symmetry with pinv and mldivide for some reason, so
    if issymmetric(A)
        iA = 0.5 * (iA + iA');
    end
end