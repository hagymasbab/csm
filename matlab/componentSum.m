function C = componentSum(coeffs,comps)
    C = zeros(size(comps{1}));
    for i=1:size(coeffs,1)
        C = C + coeffs(i,1) * comps{i};
    end
end