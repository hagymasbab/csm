function sc = cholcell(ch)
    k = size(ch,2);
    sc = cell(1,k);
    for i=1:k
        sc{i} = ch{i}' * ch{i};
    end  
end