function red = gestaltReduce(ge,new_k,new_Dv)
    red.X = ge.X;
    red.obsVar = ge.obsVar;
    red.Dx = ge.Dx;
    red.k = new_k;
    red.sparsity = ge.sparsity;
    for i=1:red.k
        temp = ge.cc{i};
        red.cc{i} = temp(1:new_Dv,1:new_Dv);
    end
    red.Dv = new_Dv;
    red.A = ge.A(:,1:red.Dv);
end
    
    