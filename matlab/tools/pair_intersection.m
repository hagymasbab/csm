function is = pair_intersection(pairs1,pairs2)
    % n x 2 matrices
    is = [];
    for p = 1:size(pairs1,1)
        f1 = find(pairs2(:,1) == pairs1(p,1));
        f2 = find(pairs2(:,2) == pairs1(p,2));
        if ~isempty(intersect(f1,f2))
            is = [is;pairs1(p,:)];
        end
    end
end