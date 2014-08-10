function chosen = chooseKfromN(k,n)
    allnum = 1:n;
    chosen = zeros(1,k);
    for i=1:k
        randidx = randi(size(allnum,2));
        chosen(1,i) = allnum(1,randidx);
        allnum(:,randidx) = [];
    end
end