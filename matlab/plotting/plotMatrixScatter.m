function plotMatrixScatter(a,b)
    if ~iscell(a)
        a = {a};
    end
    if ~iscell(b)
        a = {b};
    end
    
    a_sum = componentSum(1,a);
    b_sum = componentSum(1,b);
    
    a_all = cell2mat(a);
    b_all = cell2mat(b);
    
    