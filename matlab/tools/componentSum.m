function C = componentSum(coeffs,comps)
    if length(coeffs) == 1
        coeffs = coeffs * ones(size(comps,2));
    elseif size(coeffs,1) ~= size(comps,2)
        error('Number of coefficients (%d) is not equal to the numer of components (%d).\n',size(coeffs,1),size(comps,2));
    end
    C = zeros(size(comps{1}));
    for i=1:size(coeffs,1)
        C = C + coeffs(i,1) * comps{i};
    end
end