function [permA,permCC] = filterPermutation(A,cc,filename)
    load(filename); % orients, lambdas, maxX, maxY
    
    % cluster orientations into 4 groups
    cl_or = kmeans(orients,4);
        
    [~,sortidx] = sortrows([cl_or maxX' maxY' lambdas],[1 2 3 4]);
    permA = A(:,sortidx);
    
    permCC = {};
    for i = 1:size(cc,2)
        rowswapped = cc{i}(sortidx,:);
        permCC{i} = rowswapped(:,sortidx);
    end
end