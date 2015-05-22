function cc = eyes(k,dim,coeff)
    cc = cell(1,k);
    for i=1:k
        cc{i} = coeff * eye(dim);
%         pos = i * 8;
%         cc{i}(pos,pos+1) = 1; 
%         cc{i}(pos+1,pos) = 1;
    end
end